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
    chain_index LONG,          -- position along the chain
    position LONG,             -- bp coordinate of the second base
    fork LONG,                 -- 0 linear, else 1 or 2 in replication order
    pair SYMBOL,               -- the dinucleotide, e.g. 'AG'
    transition_index LONG,     -- 0..15
    value DOUBLE,              -- dG of the step (kcal/mol)
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
| `genome_chain.py` | the data source — FASTA reader, energy table, step events |
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

| step_index | fork | chain_index | position | pair | transition_index | value | timestamp |
|---|---|---|---|---|---|---|---|
| 1 | `fork1` | 3,925,745 | 3,925,746 | AT | 1 | -0.88 | 14:27:49.216611 |
| 2 | `fork2` | 3,925,744 | 3,925,745 | TC | 7 | -1.3 | 14:27:49.217111 |
| 3 | `fork1` | 3,925,746 | 3,925,747 | TC | 7 | -1.3 | 14:27:49.217611 |
| 4 | `fork2` | 3,925,743 | 3,925,744 | CT | 13 | -1.28 | 14:27:49.218111 |
| 5 | `fork1` | 3,925,747 | 3,925,748 | CT | 13 | -1.28 | 14:27:49.218611 |
| 6 | `fork2` | 3,925,742 | 3,925,743 | TT | 5 | -1.0 | 14:27:49.219111 |
| 7 | `fork1` | 3,925,748 | 3,925,749 | TA | 4 | -0.58 | 14:27:49.219611 |
| 8 | `fork2` | 3,925,741 | 3,925,742 | TC | 7 | -1.3 | 14:27:49.220111 |

`chain_index` moves **outward in both directions** from 3,925,744, and no two
consecutive rows are neighbours on the chain — they belong to different forks.

Filtering to one fork recovers a contiguous walk at the 1 ms replication step:

```sql
SELECT * FROM genome_steps WHERE fork = 'fork1' ORDER BY step_index;
```

| step_index | chain_index | pair | value | timestamp |
|---|---|---|---|---|
| 1 | 3,925,745 | AT | -0.88 | 14:27:49.216611 |
| 3 | 3,925,746 | TC | -1.3 | 14:27:49.217611 |
| 5 | 3,925,747 | CT | -1.28 | 14:27:49.218611 |
| 7 | 3,925,748 | TA | -0.58 | 14:27:49.219611 |
| 9 | 3,925,749 | AT | -0.88 | 14:27:49.220611 |

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
   positive ΔG (should be 0): 0
   duplicate step_index: 0
   repeated chain_index: 0
   chain gaps: 0
   fork imbalance (>1): 0
✅ Data quality PASSED
```
