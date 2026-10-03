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
