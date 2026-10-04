# SKA RealTime Genomics

**Entropy-driven, real-time learning on a genomic chain with the Structured
Knowledge Accumulation (SKA) framework.**

The chromosome is streamed one base-pair step at a time, at the rate the cell
builds it, and a forward-only learner consumes it live.

## The chain as a data stream

```
step k  =  b_k -> b_{k+1}
```

Each step carries the nearest-neighbour stacking free energy of that step
(SantaLucia 1998). Emission follows **replication order**: two forks leaving
`oriC` in opposite directions at ~1,000 bp/s each, the rate of DNA Pol III, so
one pass over E. coli K-12 MG1655 takes ~39 minutes — matching the known ~40 min
C-period.

## Every transition carries its own level

ΔG is a duplex property, and a duplex step is identical to its reverse
complement: 16 transitions, only 10 distinct values. The missing distinction is
geometric — of the six base-pair step parameters, **tilt and shift change sign**
between a step and its complement, while twist, roll, slide and rise do not.

```
level = ΔG × sign(tilt)  +  0.12 × z(twist)
```

16 distinct levels, minimum gap 0.0825. Feeding `level` instead of a return
between consecutive steps raises band separation from η² = 0.163 to **0.383**,
and each of the 16 transitions resolves into its own band of the transition
probability `P = exp(-|ΔH/H|)`.

Each band then splits in four — the previous step must end on the current step's
first letter, so there are exactly four predecessors — giving **64 sub-bands,
one per trinucleotide**. That ordering is uncorrelated with 3-mer frequency
(ρ = +0.01) and with GC content (ρ = −0.03), so it is not a re-encoding of any
standard composition statistic.

**What is not yet established:** whether band occupancy along the chromosome
tracks anything biological — coding vs intergenic, leading vs lagging strand,
the ori→ter skew inversion. That test is open.

## Layout

| folder | role |
|---|---|
| `genome-data-validation/` | the collector: FASTA → paced stream → QuestDB, with its own README |
| `transition-bands/` | the 16-band / 64-sub-band result — data, figures and how they were produced |
| `papers/` | Olson et al. 1998, Lankaš et al. 2003 |

## Quick start

```bash
cd genome-data-validation
pip install -r requirements.txt

# fetch the E. coli K-12 MG1655 reference genome
./fetch_data.sh

# stream the chain into QuestDB at the biological fork rate (~39 min)
python genome_stream_validator.py --order replication --origin 3925744 --rate 2000

# QC the collection
python validate_data.py
```

`--max-steps N` collects a subset first if you want a quick check. Ingestion
uses QuestDB's ILP on port 9009 (~600,000 rows/s); pg-wire on 8812 caps near
3,500 rows/s and cannot sustain the stream.

Each collected step carries `value` (ΔG), `twist`, `tilt` and `level` — the
16-valued input described above — so the dataset is ready for a learner without
further encoding.

## SKA framework: open science, proprietary real-time engine

The mathematical foundation and the batch implementation are public for
verification:

- [QUANTIOTA / Arxiv — Foundation Theory](https://github.com/quantiota/Arxiv)
- *Structured Knowledge Accumulation: An Autonomous Framework for Layer-Wise
  Entropy Reduction in Neural Learning* — [arXiv:2503.13942](https://arxiv.org/abs/2503.13942)
- *Structured Knowledge Accumulation: The Principle of Entropic Least Action in
  Forward-Only Neural Learning* — [arXiv:2504.03214](https://arxiv.org/abs/2504.03214)

The real-time engine extends that foundation to continuous entropy learning on a
live stream. **That part is proprietary and is not included in this
repository** — the results in `transition-bands/` were produced with it, and the
data they were computed from is reproducible here in full.

## References

> Olson, Gorin, Lu, Hock & Zhurkin (1998), *DNA sequence-dependent deformability
> deduced from protein–DNA crystal complexes*, PNAS 95:11163–11168.
>
> Lankaš, Šponer, Langowski & Cheatham (2003), *DNA basepair step deformability
> inferred from molecular dynamics simulations*, Biophys. J. 85:2872–2883.
>
> SantaLucia (1998), *A unified view of polymer, dumbbell, and oligonucleotide
> DNA nearest-neighbor thermodynamics*, PNAS 95:1460–1465.

## Citation

> Bouarfa Mahi, *SKA RealTime Genomics: entropy-driven real-time learning on a
> genomic chain* (2026), GitHub.

## License

MIT





