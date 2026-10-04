# The 16 transition bands, and the 64 trinucleotide sub-bands

SKA run on the E. coli K-12 MG1655 replication chain, with each base-pair step
encoded so that **every one of the 16 transitions carries its own input level**.
The transition probability

```
P(k) = exp( -| (H_k - H_{k-1}) / H_k | )
```

then separates into 16 bands, each splitting into exactly 4 sub-bands — 64 in
total, one per trinucleotide.

## The run

| | |
|---|---|
| record | `NC_000913.3` — E. coli K-12 MG1655 |
| series | `fork1`, replication order from `oriC` (3,925,744) |
| steps | 3,500 — `chain_index` 3,925,745–3,929,244 |
| input | `level` mode: the step's own signed level, **no differencing** |
| scale | 1.0 |
| buffer | 3,500 — the matrix grows the whole way, nothing pruned |
| x range | 0.0938–0.8555, nothing clipped, nothing on the rails |
| entropy | −1.5750 … −0.0000, 100 % negative |

## Why each transition owns a band

ΔG is a duplex property, and a duplex step is identical to its reverse
complement, so 16 steps carry only 10 ΔG values. The missing distinction is
geometric: of the six base-pair step parameters, **tilt and shift change sign**
between a step and its complement (Olson et al. 1998 and Lankaš et al. 2003 both
state it in the notes to their Table 1). Twist breaks the remaining −1.44/−1.45
near-tie.

```
level = ΔG × sign(tilt)  +  0.12 × z(twist)      ->  16 distinct, min gap 0.0825
```

Feeding `level` directly — rather than a return between consecutive steps —
makes the input a function of the current transition alone. Measured band
separation:

| input | η² |
|---|---|
| return mode, best scale (10) | 0.163 |
| **level mode, scale 1** | **0.383** |

η² is the fraction of the variance in P explained by which transition the step
is: 0 = bands fully overlapped, 1 = fully separated.

## Why each band splits in four

Consecutive dinucleotides overlap by one letter. If the current step is `X→Y`,
the previous step must have *ended* on X, so it is one of `AX`, `TX`, `GX`, `CX`
— four predecessors, never more. P depends on the change between the two steps,
so each transition splits into exactly four sub-bands:

```
16 transitions x 4 predecessors = 64 = the trinucleotides
```

Reading `prev→cur` as `Z X Y`, the four groups under `A→A` are `AAA`, `GAA`,
`TAA`, `CAA`.

The predecessor often matters more than the transition itself:

```
G->C    CG-> 0.945    AG-> 0.117    TG-> 0.093    GG-> 0.086
T->G    TT-> 0.589    AT-> 0.573    CT-> 0.547    GT-> 0.533
```

## The ordering is not a composition statistic

Every standard 64-number table in genomics — trinucleotide frequencies, order-2
Markov conditionals, codon usage — is built by counting occurrences. This one is
not, and the correlations are flat:

```
Spearman( mean P , 3-mer count in window )  = +0.07
Spearman( mean P , 3-mer count in genome )  = +0.01
Spearman( mean P , GC content of the 3-mer) = -0.03
```

So the band positions are an axis orthogonal to composition, fixed by energy,
strand direction and twist. Because they are fixed, differences between windows
of the chain come entirely from the sequence — a fixed ruler, a varying
measurement.

What this does **not** yet establish: whether band occupancy along the
chromosome tracks anything biological (coding vs intergenic, leading vs lagging
strand, the ori→ter skew inversion). That test is open.

## Files

| file | contents |
|---|---|
| `ska_results_fork1_level_scale1.csv` | the full run, 3,500 rows, exported from QuestDB |
| `band_means.csv` | mean, min and max P for the 16 transitions and the 64 trinucleotides |
| `images/` | the panels |

### `ska_results_fork1_level_scale1.csv`

One row per step. Note two columns whose names come from the engine's return
mode and read differently here:

| column | meaning in level mode |
|---|---|
| `value` | the step's ΔG, kcal/mol — carried through, not the input |
| `value_return` | **the `level`**, i.e. the actual feature fed to the sigmoid |
| `x_input` | `sigmoid(level × scale)` — what the learner received |
| `entropy` | H(k) |
| `knowledge`, `decision`, `decision_norm` | ‖Z‖, D[-1], ‖D‖ |
| `matrix_size` | k, growing to 3,500 |

The export has no `pair` column; `band_means.csv` carries the per-transition
figures, and the Grafana queries join `genome_steps` for the letters.

## Images to add

| file | panel |
|---|---|
| `images/entropy.png` | entropy vs `step_index`, all 3,500 steps |
| `images/probability_16_bands.png` | P vs `step_index`, partitioned by transition |
| `images/band_AA.png` … `images/band_CC.png` | one per transition, partitioned by trinucleotide — 16 files |

Naming for the 16: `band_<pair>.png`, e.g. `band_AA.png`, `band_GC.png`.

Queries that produce them are in `../transition-bands/queries/`:

- `grafana_entropy.sql` — the entropy panel
- `grafana_transition_bands.sql` — the 16 bands


**Fix the axes before capturing the set of 16** — y min 0, y max 1, x max 7000.
Otherwise each panel autoscales and unequal spreads look equal: `G→C` spans
0.09–0.95 while `T→G` spans only 0.53–0.59.

## References

> Olson, Gorin, Lu, Hock & Zhurkin (1998), *DNA sequence-dependent deformability
> deduced from protein–DNA crystal complexes*, PNAS 95:11163–11168 — Table 1.
>
> Lankaš, Šponer, Langowski & Cheatham (2003), *DNA basepair step deformability
> inferred from molecular dynamics simulations*, Biophys. J. 85:2872–2883.
