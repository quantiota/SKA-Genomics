"""Genomic chain reader — the data source for the validation stream.

Turns a FASTA record into a sequence of samples, where a sample is one
base-pair step along the chain.

    step k  =  b_k -> b_{k+1}
    value   =  dG37 of that step   (SantaLucia 1998 nearest-neighbour table)

The value sits on the STEP, not on the letter: a base is a token with no
magnitude, a step is a measured physical quantity.
"""
from __future__ import annotations

import numpy as np

BASES = "ATGC"                              # A=0, T=1, G=2, C=3

# dG37 (kcal/mol), SantaLucia 1998 unified nearest-neighbour parameters,
# indexed by 4*first + next. Sixteen steps, ten distinct values: complementary
# steps (AA/TT, CA/TG, GT/AC, CT/AG, GA/TC, GG/CC) share a parameter because
# they are the same physical stack read from opposite strands.
ENERGY = np.array([
    -1.00, -0.88, -1.28, -1.44,             # A->A  A->T  A->G  A->C
    -0.58, -1.00, -1.45, -1.30,             # T->A  T->T  T->G  T->C
    -1.30, -1.44, -1.84, -2.24,             # G->A  G->T  G->G  G->C
    -1.45, -1.28, -2.17, -1.84,             # C->A  C->T  C->G  C->C
])

# ---- the geometry that breaks the ten-value degeneracy ----------------
#
# A duplex step is identical to its reverse complement -- the dyad rotation
# maps one onto the other -- so dG gives ten values for sixteen steps. Of the
# six base pair step parameters, four (twist, roll, slide, rise) are invariant
# under that rotation and two (TILT and SHIFT) change sign. Both source papers
# state it in the notes to their Table 1:
#
#   Olson et al. 1998, PNAS 95:11163 -- "AA and TT, AG and CT, etc., have
#   identical averages except for different signs of Tilt and Shift."
#   Lankas et al. 2003, Biophys J 85:2872 -- tilt and shift "change sign upon
#   changing the direction in which a DNA sequence is followed."
#
# sign(tilt) is therefore a MEASURED quantity, not a convention: it is the one
# bit that says which strand the step is read from. Twist is added at a small
# weight to break the dG near-tie between AC/GT (-1.44) and CA/TG (-1.45),
# whose twists are 31.5 and 37.3 deg.
#
#     LEVEL = dG * sign(tilt)  +  TWIST_WEIGHT * z(twist)
#
# TWIST_WEIGHT = 0.12 maximises the worst-case separation: sixteen distinct
# levels, minimum gap 0.0825 against 0.010 for signed dG alone. It is sharply
# tuned -- 0.14 drops the gap to 0.03 as accidental collisions reappear.
#
# Tilt magnitudes (0.1-1.7 deg) sit far below their dispersion (~3 deg), so
# only the SIGN is used. For GG/CC and AC/GT it rests on +-0.1 deg and is a
# weak mean tendency; for AA, AG and GA (1.4-1.7 deg) it is firm.
#
# Olson et al. 1998, Table 1 (protein-DNA crystal complexes), in the same
# 4*first + next order as ENERGY. Twist is strand-symmetric; tilt is
# antisymmetric, so a step and its reverse complement carry opposite signs and
# the four self-complementary steps (AT, TA, GC, CG) are exactly zero.
TWIST = np.array([
     35.1,  29.3,  31.9,  31.5,             # A->A  A->T  A->G  A->C
     37.8,  35.1,  37.3,  36.3,             # T->A  T->T  T->G  T->C
     36.3,  31.5,  32.9,  33.6,             # G->A  G->T  G->G  G->C
     37.3,  31.9,  36.1,  32.9,             # C->A  C->T  C->G  C->C
])
TILT = np.array([
     -1.4,   0.0,  -1.7,  -0.1,             # A->A  A->T  A->G  A->C
      0.0,  +1.4,  -0.5,  +1.5,             # T->A  T->T  T->G  T->C
     -1.5,  +0.1,  -0.1,   0.0,             # G->A  G->T  G->G  G->C
     +0.5,  +1.7,   0.0,  +0.1,             # C->A  C->T  C->G  C->C
])
TWIST_WEIGHT = 0.12

LEVEL = (ENERGY * np.where(TILT < 0, -1.0, 1.0)
         + TWIST_WEIGHT * (TWIST - TWIST.mean()) / TWIST.std())


def read_fasta(path):
    """Return (header, sequence) of the first record in a FASTA file."""
    header, chunks = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    break
                header = line[1:]
            elif line:
                chunks.append(line.upper())
    return header, "".join(chunks)


# reverse complement of each step index, in the same 4*first+next encoding.
# A->T=1, T->A=0, G->C=3, C->G=2 under BASES = "ATGC".
_COMP = {0: 1, 1: 0, 2: 3, 3: 2}
RC_INDEX = np.array([4 * _COMP[i % 4] + _COMP[i // 4] for i in range(16)], dtype=np.int16)
RC_PAIR = {a + b: BASES[_COMP[BASES.index(b)]] + BASES[_COMP[BASES.index(a)]]
           for a in BASES for b in BASES}


def transition_index(seq: str, circular: bool = False) -> np.ndarray:
    """Index 0-15 for every step i -> i+1; -1 where a base is not ACGT.

    circular: also emit the closing step last base -> first base, so a circular
    chromosome of N bases yields N steps rather than N-1.
    """
    lut = np.full(256, -1, dtype=np.int16)
    for code, b in enumerate(BASES):
        lut[ord(b)] = code
    c = lut[np.frombuffer(seq.encode(), dtype=np.uint8)]
    first = np.concatenate([c, c[:1]])[:-1] if circular else c[:-1]
    nxt = np.concatenate([c[1:], c[:1]]) if circular else c[1:]
    ok = (first >= 0) & (nxt >= 0)
    idx = np.full(first.size, -1, dtype=np.int16)
    idx[ok] = 4 * first[ok] + nxt[ok]
    return idx


class GenomeChain:
    """Walks a FASTA record and yields one event per base-pair step."""

    def __init__(self, fasta_path, max_steps=None, order="linear", origin=None,
                 circular=None, leading_strand=True):
        """
        order           'linear'      position 1 -> end
                        'replication' two forks outward from oriC
        circular        close the ring (last base -> first base). Defaults to
                        True for replication order, where the forks must wrap.
        leading_strand  in replication order, report fork 2's steps as the
                        REVERSE COMPLEMENT, i.e. as its own leading strand reads
                        them. dG is unaffected (it is strand-symmetric); only
                        `pair` and `transition_index` change. With this on, GC
                        skew keeps one sign across the whole stream instead of
                        reversing at oriC and ter.
        """
        if order not in ("linear", "replication"):
            raise ValueError("order must be 'linear' or 'replication'")
        self.header, self.seq = read_fasta(fasta_path)
        self.circular = (order == "replication") if circular is None else bool(circular)
        self.idx = transition_index(self.seq, circular=self.circular)
        self.max_steps = max_steps
        self.order = order
        self.origin = origin
        self.leading_strand = bool(leading_strand)

    @property
    def accession(self) -> str:
        """Record ID taken from the FASTA header, not assumed."""
        return self.header.split()[0] if self.header else ""

    @property
    def chain_length(self) -> int:
        return len(self.seq)

    @property
    def total_steps(self) -> int:
        return int(self.idx.size)

    def stats(self) -> dict:
        """Whole-chain statistics (evaluation only, never model input)."""
        g = ENERGY[self.idx[self.idx >= 0]]
        return {
            "header": self.header,
            "chain_length": self.chain_length,
            "total_steps": self.total_steps,
            "value_mean": float(g.mean()),
            "value_std": float(g.std()),
            "value_min": float(g.min()),
            "value_max": float(g.max()),
            "value_total": float(g.sum()),
        }

    def _emit(self, k, seq_no, fork=0):
        i = int(self.idx[k])
        if i < 0:
            return None                       # skip steps touching non-ACGT
        N = len(self.seq)
        pair = self.seq[k] + self.seq[(k + 1) % N]
        if fork == 2 and self.leading_strand:
            i = int(RC_INDEX[i])              # fork 2 reads the other strand
            pair = RC_PAIR[pair]
        return {
            "step": seq_no,                   # emission order
            "chain_index": k + 1,             # position along the chain
            "position": (k + 1) % N + 1,      # bp coordinate of the second base
            "pair": pair,
            "index": i,
            "value": float(ENERGY[i]),        # identical either strand
            "twist": float(TWIST[i]),         # identical either strand
            "tilt": float(TILT[i]),           # OPPOSITE on the other strand
            "level": float(LEVEL[i]),         # 16 distinct, one per transition
            "fork": fork,
        }

    def _order_indices(self):
        """Chain indices in emission order.

        linear       0, 1, 2, ... — a plain walk, NOT the order the cell reads
        replication  two forks leaving oriC in opposite directions around the
                     circular chromosome, interleaved, meeting at the terminus
        """
        n_lim = self.total_steps if self.max_steps is None \
            else min(self.max_steps, self.total_steps)
        if self.order == "linear":
            yield from ((k, 0) for k in range(n_lim))
            return
        N = self.total_steps
        o = (self.origin if self.origin is not None else 0) % N
        emitted = 0
        for d in range(1, N // 2 + 2):
            for fork, k in ((1, (o + d - 1) % N), (2, (o - d) % N)):
                if emitted >= n_lim:
                    return
                yield k, fork
                emitted += 1

    def events(self):
        """Yield one event per step, in emission order."""
        seq_no = 0
        for k, fork in self._order_indices():
            seq_no += 1
            e = self._emit(k, seq_no, fork)
            if e is not None:
                yield e


if __name__ == "__main__":
    import sys
    chain = GenomeChain(sys.argv[1])
    for k, v in chain.stats().items():
        print(f"{k:14s} {v}")
    print("\nfirst 5 steps:")
    for e in list(chain.events())[:5]:
        print("  ", e)
