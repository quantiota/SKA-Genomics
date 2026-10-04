# Genomic Data Stream Validation

Simulated real-time genomic chain streaming and QuestDB ingestion validation.

## Purpose

- Validate streaming of a genomic chain one base-pair step at a time
- Test QuestDB table creation and ILP ingestion
- Collect clean raw genomic parameters for future SKA analysis

No returns, no entropy, no learning at this stage. **Raw values only.**

## Development environment

AI Agent Host — [Quick start](https://github.com/quantiota/AI-Agent-Host)

## Quick start

```bash
pip install -r requirements.txt

# Download the E. coli K-12 MG1655 reference genome from NCBI
./fetch_data.sh

# Point at your QuestDB instance
export QDB_PG_HOST=localhost      # pg-wire on 8812, ILP on 9009

# 1. Test the connection
python test_connection.py

# 2. Stream the chain into QuestDB
python genome_stream_validator.py

# 3. Validate the collected data
python validate_data.py
```

Options: `--max-steps N` for a subset, `--rate HZ` for the total emission rate,
`--order linear|replication`, `--origin N` for oriC, `--raw-strand` for
top-strand labels on fork 2, and `--recreate` to drop an existing non-empty
table (refused without it, so a running collection cannot be destroyed by
accident).

Note: interrupting during a socket send may leave a partially written packet;
`validate_data.py` will report the resulting gap.

## What a step is

The chain is walked one base-pair step at a time. Each step carries the
**nearest-neighbour stacking free energy** of that step:

```
step k  =  b_k -> b_{k+1}
value   =  dG37 of that step     (SantaLucia 1998)
```

Sixteen dinucleotide steps, ten distinct values — complementary steps
(`AA/TT`, `CA/TG`, `GT/AC`, `CT/AG`, `GA/TC`, `GG/CC`) share a parameter,
because they are the same physical stack read from opposite strands.

The value sits on the **step**, not on the letter: a base is a token with no
magnitude, a step is a measured physical quantity.

### The ten-value degeneracy, and the geometry that breaks it

A duplex step is **the same physical object** as its reverse complement: the
dyad rotation that swaps the two strands maps `AA` onto `TT`. Sixteen names,
ten objects. ΔG cannot separate a pair — not because the measurement is coarse,
but because there is nothing there to separate.

Of the six base-pair step parameters, four are invariant under that rotation
(twist, roll, slide, rise) and **two change sign** — tilt and shift. Both source
papers say so in the notes to their Table 1:

> "AA and TT, AG and CT, etc., have identical averages **except for different
> signs of Tilt and Shift**." — Olson et al. 1998
>
> tilt and shift "**change sign** upon changing the direction in which a DNA
> sequence is followed" — Lankaš et al. 2003

So `sign(tilt)` is the one bit that says which strand a step is read from, and
it is **measured**, not conventional. Twist is added at a small weight to break
the ΔG near-tie between `AC`/`GT` (−1.44) and `CA`/`TG` (−1.45), whose twists
are 31.5° and 37.3°:

```
level = ΔG × sign(tilt)  +  0.12 × z(twist)
```

That gives **sixteen distinct levels, one per transition**, minimum gap 0.0825
against 0.010 for signed ΔG alone. The weight 0.12 maximises the worst-case
separation and is sharply tuned — 0.14 drops the gap to 0.03 as accidental
collisions reappear.

Each step therefore carries `twist`, `tilt` and `level` alongside `value`.
Note `tilt` is the one column that is **not** strand-invariant: fork 2's steps
carry the opposite sign, which is the physically correct answer and the reason
the degeneracy lifts.

Caveats worth keeping in view: tilt magnitudes (0.1–1.7°) sit far below their
dispersion (~3°), so only the **sign** is used, and for `GG`/`CC` and `AC`/`GT`
that sign rests on ±0.1° — a weak mean tendency. For `AA`, `AG` and `GA`
(1.4–1.7°) it is firm.

> Olson, Gorin, Lu, Hock & Zhurkin (1998), *DNA sequence-dependent deformability
> deduced from protein–DNA crystal complexes*, PNAS 95:11163–11168 — Table 1.
> Lankaš, Šponer, Langowski & Cheatham (2003), *DNA basepair step deformability
> inferred from molecular dynamics simulations*, Biophys. J. 85:2872–2883.

## Why 1,000 steps per second

A genomic chain has no native sampling rate, so one must be imposed. Rather than
pick a convenient number, the stream runs at the rate the cell builds the chain:

```
E. coli replication fork (DNA Pol III)  ~1,000 bp/s   ->  1 ms per step
```

Cross-check: replication is bidirectional from `oriC`, so each of two forks
covers ~2.32 Mb, and the known ~40 min C-period gives
`2.32e6 / 2400 s ~ 967 bp/s`.

One single pass over 4,641,651 steps therefore takes **~77 minutes**.

## Timestamps are assigned, not measured

Each row is stamped on a **fixed grid** anchored at the run's start time:

```
timestamp(k) = t_start + k / rate_hz
```

These are **assigned** timestamps. A downstream `delta_t` will return exactly
`1/rate_hz` with no jitter — it recovers the configured rate rather than
measuring anything. The determinism is deliberate and is the right property for
analysis, but it must not be mistaken for an observation.

Packets are paced against the wall clock so the collection takes real time, but
rows are shipped 1,000 at a time, not at the instant each timestamp claims. To
obtain genuine jitter the stream would have to stamp every row with
`time.time_ns()` at the moment of emission, one row at a time.

In replication order `--rate` is the total, shared between two forks, so the
grid spacing differs depending on how you read the table:

```
raw stream, forks interleaved   delta_t = 0.5 ms  (at --rate 2000)
one fork alone                  delta_t = 1.0 ms  = 1,000 bp/s
```

The 0.5 ms is an artefact of two forks sharing one emission channel; nothing
physical moves at that rate. Filtering by `fork` gives the biological 1 ms.

## Emission order

```
linear       position 1 -> end                                    (default)
replication  two forks leaving oriC in opposite directions
             around the circular chromosome, interleaved
```

**`linear` is not how the cell reads the chain.** Replication starts at `oriC`
(~3.92 Mb in E. coli) and runs as two forks toward the terminus. This matters:
the GC-skew dipole reverses at `oriC` and at `ter`, so a linear walk crosses
both reversals at positions that look arbitrary to a learner.

```bash
python genome_stream_validator.py --order replication --origin 3925744 --rate 2000
```

`--origin` is a **0-based** chain index, so fork 1 begins at `chain_index`
3,925,745 and fork 2 at 3,925,744. Replication order closes the ring — the
step from the last base back to the first is emitted — giving 4,641,652 steps
against the linear 4,641,651.

### Read one fork at a time

In replication order the two forks are **interleaved**, so consecutive rows sit
about 2 Mb apart on the chain and **no two consecutive rows are neighbours**.
That is physically correct — both forks advance simultaneously — but a consumer
that reads transitions from consecutive rows would be reading noise.

`fork` is a SYMBOL column for exactly this reason. Read one series at a time:

```sql
SELECT * FROM genome_steps WHERE fork = 'fork1' ORDER BY step_index;
```

The same applies to timing. At `--rate 2000` consecutive raw rows are 0.5 ms
apart, but consecutive rows *within* a fork are 1.0 ms apart — the replication
step. Reading the raw interleave gets both the positions and the intervals
wrong.

### The rate is shared between the forks

`--rate` is the **total** emission rate. In replication order each fork
therefore advances at half of it:

| `--rate` | per fork | duration |
|---|---|---|
| 1000 | 500 bp/s | 77 min |
| **2000** | **1,000 bp/s** | **39 min** |

`--rate 2000` gives each fork the biological ~1,000 bp/s, and the resulting
~39 min matches the known ~40 min C-period — a useful consistency check.

### Which strand fork 2 is reported on

Fork 2 travels backwards along top-strand coordinates, so its own leading strand
is the reverse complement. By default each step is reported **as its own fork
reads it**, which keeps GC skew at one sign across the whole stream rather than
reversing at `oriC` and `ter`. Pass `--raw-strand` for raw top-strand labels.

This affects `pair` and `transition_index` only. **`value` is identical either
way**, because dG is strand-symmetric by construction — `AA/TT`, `CA/TG`,
`GT/AC`, `CT/AG`, `GA/TC` and `GG/CC` each share a parameter. The choice is
about metadata, not about what a learner consumes.

## Schema

```sql
CREATE TABLE genome_steps (
    record_id SYMBOL,          -- NC_000913.3
    organism SYMBOL,
    step_index LONG,           -- emission order
    fork SYMBOL,               -- 'linear', or 'fork1' / 'fork2' in replication order
    chain_index LONG,          -- position along the chain
    position LONG,             -- bp coordinate of the second base
    pair SYMBOL,               -- the dinucleotide, e.g. 'AG'
    transition_index LONG,     -- 0..15
    value DOUBLE,              -- dG of the step (kcal/mol), strand-symmetric
    twist DOUBLE,              -- deg, Olson 1998, strand-symmetric
    tilt DOUBLE,               -- deg, Olson 1998, OPPOSITE on the other strand
    level DOUBLE,              -- dG*sign(tilt) + 0.12*z(twist): 16 distinct
    rate_hz DOUBLE,            -- imposed acquisition rate (provenance)
    chain_length LONG,
    total_steps LONG,
    timestamp TIMESTAMP        -- assigned, on a fixed grid
) TIMESTAMP(timestamp) PARTITION BY DAY;
```

Ingestion uses QuestDB's **InfluxDB line protocol on port 9009**
(~600,000 rows/s), not pg-wire — which caps near 3,500 rows/s and cannot
sustain the stream.

## Files

| file | role |
|---|---|
| `genome_chain.py` | the data source — FASTA reader, energy and geometry tables, step events |
| `genome_stream_validator.py` | paced emission into QuestDB over ILP |
| `test_connection.py` | QuestDB smoke test |
| `validate_data.py` | QC report on the collected chain |
| `config.py` | QuestDB, ILP, genome and logging configuration |
| `fetch_data.sh` | download the reference genome from NCBI |

## Output

- `logs/genome_validation.log`

## Data Collection Summary

Produced by `python validate_data.py` after a full replication-order collection.

```
record     NC_000913.3  Escherichia coli K-12 MG1655
chain      4,641,652 bp (circular)
collected  4,641,652 steps  (order=replication, oriC 3,925,744,
           2,000 Hz total = 1,000 bp/s per fork, 38.7 min)
```

The 38.7 min is worth noting against the known ~40 min C-period: the rate was
set from the fork speed, not fitted to the duration.

### Chain coverage

| | |
|---|---|
| rows | 4,641,652 |
| distinct `chain_index` | 4,641,652 |
| range | 1–4,641,652 |
| repeated positions | 0 |
| gaps within range | 0 |

The step count is 4,641,652 rather than 4,641,651 because replication order
closes the ring: the step from the last base back to the first is emitted.

### Fork balance

| fork | steps |
|---|---|
| `fork1` | 2,320,826 |
| `fork2` | 2,320,826 |
| difference | 0 |

Each fork covers exactly half the chromosome, which is the check that the
interleaving is correct.

### Collected rows

The first eight rows as collected — the two forks interleaved, leaving `oriC`
in opposite directions:

| step_index | fork | chain_index | position | pair | transition_index | value | tilt | level | timestamp |
|---|---|---|---|---|---|---|---|---|---|
| 1 | `fork1` | 3,925,745 | 3,925,746 | AT | 1 | -0.88 | +0.0 | -1.1150 | 10:41:40.935027 |
| 2 | `fork2` | 3,925,744 | 3,925,745 | TC | 7 | -1.3 | +1.5 | -1.1976 | 10:41:40.935527 |
| 3 | `fork1` | 3,925,746 | 3,925,747 | TC | 7 | -1.3 | +1.5 | -1.1976 | 10:41:40.936027 |
| 4 | `fork2` | 3,925,743 | 3,925,744 | CT | 13 | -1.28 | +1.7 | -1.3897 | 10:41:40.936527 |
| 5 | `fork1` | 3,925,747 | 3,925,748 | CT | 13 | -1.28 | +1.7 | -1.3897 | 10:41:40.937027 |
| 6 | `fork2` | 3,925,742 | 3,925,743 | TT | 5 | -1.0 | +1.4 | -0.9554 | 10:41:40.937527 |
| 7 | `fork1` | 3,925,748 | 3,925,749 | TA | 4 | -0.58 | +0.0 | -0.4052 | 10:41:40.938027 |
| 8 | `fork2` | 3,925,741 | 3,925,742 | TC | 7 | -1.3 | +1.5 | -1.1976 | 10:41:40.938527 |

`chain_index` moves **outward in both directions** from 3,925,744, and no two
consecutive rows are neighbours on the chain — they belong to different forks.

Filtering to one fork recovers a contiguous walk at the 1 ms replication step:

```sql
SELECT * FROM genome_steps WHERE fork = 'fork1' ORDER BY step_index;
```

| step_index | chain_index | pair | value | tilt | level | timestamp |
|---|---|---|---|---|---|---|
| 1 | 3,925,745 | AT | -0.88 | +0.0 | -1.1150 | 10:41:40.935027 |
| 3 | 3,925,746 | TC | -1.3 | +1.5 | -1.1976 | 10:41:40.936027 |
| 5 | 3,925,747 | CT | -1.28 | +1.7 | -1.3897 | 10:41:40.937027 |
| 7 | 3,925,748 | TA | -0.58 | +0.0 | -0.4052 | 10:41:40.938027 |
| 9 | 3,925,749 | AT | -0.88 | +0.0 | -1.1150 | 10:41:40.939027 |

`step_index` advances by 2 (the other fork takes the alternate slots),
`chain_index` by 1, and the timestamps by exactly 1.000 ms.

### Transition coverage — 16 of 16 present

| idx | step | ΔG | count |
|---|---|---|---|
| 0 | AA | -1.00 | 340,209 |
| 1 | AT | -0.88 | 309,950 |
| 2 | AG | -1.28 | 240,989 |
| 3 | AC | -1.44 | 246,218 |
| 4 | TA | -0.58 | 212,024 |
| 5 | TT | -1.00 | 337,381 |
| 6 | TG | -1.45 | 338,864 |
| 7 | TC | -1.30 | 258,490 |
| 8 | GA | -1.30 | 276,289 |
| 9 | GT | -1.44 | 266,254 |
| 10 | GG | -1.84 | 289,428 |
| 11 | GC | -2.24 | 384,102 |
| 12 | CA | -1.45 | 308,843 |
| 13 | CT | -1.28 | 233,173 |
| 14 | CG | -2.17 | 346,793 |
| 15 | CC | -1.84 | 252,645 |

Labels are leading-strand: fork 2's steps are reverse-complemented, so the
counts are not the raw top-strand counts. `value` is unaffected.

### Step geometry — level separation

Each transition and the level it carries, ordered by `level`. The sign is
`sign(tilt)`; the small offsets within a sign group come from the twist term.

| step | ΔG | twist ° | tilt ° | level | count |
|---|---|---|---|---|---|
| G→C | -2.24 | 33.6 | +0.0 | **-2.2677** | 384,102 |
| C→G | -2.17 | 36.1 | +0.0 | **-2.0772** | 346,793 |
| C→C | -1.84 | 32.9 | +0.1 | **-1.9015** | 252,645 |
| G→T | -1.44 | 31.5 | +0.1 | **-1.5690** | 266,254 |
| C→T | -1.28 | 31.9 | +1.7 | **-1.3897** | 233,173 |
| C→A | -1.45 | 37.3 | +0.5 | **-1.2993** | 308,843 |
| T→C | -1.30 | 36.3 | +1.5 | **-1.1976** | 258,490 |
| A→T | -0.88 | 29.3 | +0.0 | **-1.1150** | 309,950 |
| T→T | -1.00 | 35.1 | +1.4 | **-0.9554** | 337,381 |
| T→A | -0.58 | 37.8 | +0.0 | **-0.4052** | 212,024 |
| A→A | -1.00 | 35.1 | -1.4 | **+1.0446** | 340,209 |
| A→G | -1.28 | 31.9 | -1.7 | **+1.1703** | 240,989 |
| A→C | -1.44 | 31.5 | -0.1 | **+1.3110** | 246,218 |
| G→A | -1.30 | 36.3 | -1.5 | **+1.4024** | 276,289 |
| T→G | -1.45 | 37.3 | -0.5 | **+1.6007** | 338,864 |
| G→G | -1.84 | 32.9 | -0.1 | **+1.7785** | 289,428 |

```
16 distinct levels for 16 transitions   min gap 0.0825
✅ every transition has its own level
   self-complementary tilt = 0 (AT, TA, GC, CG): ✅
   tilt antisymmetric across the 6 complementary pairs: ✅
```

The six ΔG collisions are gone: `A→A` and `T→T` both carry −1.00 but sit at
+1.0446 and −0.9554, and the −1.44/−1.45 near-tie (`A→C`/`T→G`) is separated by
twist to +1.3110 and +1.6007.

### Collected vs whole-genome counts

An independent round-trip: counts read back from the database against the
reference re-read from the FASTA. All ten energy levels agree.

| ΔG | collected | genome |
|---|---|---|
| -2.24 | 384,102 | 384,102 |
| -2.17 | 346,793 | 346,793 |
| -1.84 | 542,073 | 542,073 |
| -1.45 | 647,707 | 647,707 |
| -1.44 | 512,472 | 512,472 |
| -1.30 | 534,779 | 534,779 |
| -1.28 | 474,162 | 474,162 |
| -1.00 | 677,590 | 677,590 |
| -0.88 | 309,950 | 309,950 |
| -0.58 | 212,024 | 212,024 |

### Data quality

```
🔍 Data quality:
   NULL value: 0
   NULL record_id: 0
   NULL step_index: 0
   NULL pair: 0
   NULL twist: 0
   NULL tilt: 0
   NULL level: 0
   positive ΔG (should be 0): 0
   duplicate step_index: 0
   repeated chain_index: 0
   chain gaps: 0
   fork imbalance (>1): 0
✅ Data quality PASSED
```
