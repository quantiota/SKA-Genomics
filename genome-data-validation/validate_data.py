"""Validate collected genomic data in QuestDB."""
import collections

import psycopg2

from config import QDB_CONFIG, GENOME_CONFIG
from genome_chain import ENERGY, GenomeChain

TABLE = "genome_steps"


def validate_genome_data():
    """Run validation queries on collected genomic step data."""
    try:
        conn = psycopg2.connect(**QDB_CONFIG)
        with conn.cursor() as cur:
            cur.execute("SHOW TABLES;")
            if TABLE not in [t[0] for t in cur.fetchall()]:
                print(f"❌ {TABLE} table not found")
                return

            cur.execute(f"SELECT COUNT(*) FROM {TABLE};")
            total = cur.fetchone()[0]
            print(f"📊 Total steps collected: {total:,}")

            cur.execute(f"""SELECT record_id, organism, COUNT(*) AS n,
                                   MIN(step_index), MAX(step_index), MAX(total_steps)
                            FROM {TABLE} GROUP BY record_id, organism;""")
            print("\n📍 Steps per record:")
            for rec, org, n, lo, hi, tot in cur.fetchall():
                pct = 100.0 * n / tot if tot else 0.0
                print(f"   {rec} ({org}): {n:,} steps  index {lo:,}–{hi:,}  "
                      f"{pct:.1f}% of {tot:,}")

            cur.execute(f"""SELECT record_id, MIN(value), MAX(value), AVG(value), STDDEV(value)
                            FROM {TABLE} GROUP BY record_id;""")
            print("\n📈 Value statistics (ΔG, kcal/mol):")
            for rec, lo, hi, avg, sd in cur.fetchall():
                sd = sd if sd is not None else float('nan')
                print(f"   {rec}: {lo:.2f} to {hi:.2f}  (avg {avg:.4f}, std {sd:.4f})")

            cur.execute(f"""SELECT transition_index, pair, COUNT(*) AS n, MIN(value)
                            FROM {TABLE} GROUP BY transition_index, pair
                            ORDER BY transition_index;""")
            trans = cur.fetchall()
            print(f"\n🔤 Transition coverage: {len(trans)} of 16 present")
            for ti, pair, n, v in trans:
                print(f"   {ti:2d}  {pair}  ΔG={v:+.2f}  count={n:,}")

            cur.execute(f"""SELECT record_id, step_index, pair, value, timestamp
                            FROM {TABLE} ORDER BY timestamp DESC LIMIT 10;""")
            print("\n🕐 Recent steps:")
            for rec, si, pair, val, ts in cur.fetchall():
                print(f"   {rec} step {si:,} {pair} ΔG={val:+.2f} at {ts}")

            # ---- chain coverage: every chain_index exactly once -----------
            cur.execute(f"""SELECT count(), count_distinct(chain_index),
                                   min(chain_index), max(chain_index)
                            FROM {TABLE};""")
            n, n_distinct, lo, hi = cur.fetchone()
            print("\n🧬 Chain coverage:")
            print(f"   rows {n:,}   distinct chain_index {n_distinct:,}   "
                  f"range {lo:,}–{hi:,}")
            print(f"   repeated positions: {n - n_distinct:,}")
            print(f"   gaps within range:  {(hi - lo + 1) - n_distinct:,}")

            # ---- fork balance ---------------------------------------------
            cur.execute(f"SELECT fork, count() FROM {TABLE} GROUP BY fork;")
            forks = dict(cur.fetchall())
            print("\n🍴 Fork balance:")
            for f, c in sorted(forks.items()):
                print(f"   {f}: {c:,}")
            fork_imbalance = 0
            if "fork1" in forks and "fork2" in forks:
                fork_imbalance = abs(forks["fork1"] - forks["fork2"])
                print(f"   |fork1 - fork2| = {fork_imbalance}  "
                      f"({'OK' if fork_imbalance <= 1 else 'UNBALANCED'})")

            # ---- transition counts vs the whole genome --------------------
            print("\n⚖️  Collected vs whole-genome counts:")
            try:
                chain = GenomeChain(GENOME_CONFIG["fasta"],
                                    circular=("fork1" in forks))
                expected = collections.Counter(
                    float(ENERGY[i]) for i in chain.idx[chain.idx >= 0])
                cur.execute(f"SELECT value, count() FROM {TABLE} GROUP BY value;")
                got = {round(v, 2): c for v, c in cur.fetchall()}
                complete = n == chain.total_steps
                if not complete:
                    print(f"   partial collection ({n:,} of {chain.total_steps:,}) "
                          f"— counts not expected to match")
                mismatches = 0
                for dG in sorted(expected):
                    e, g = expected[dG], got.get(round(dG, 2), 0)
                    flag = "" if (e == g or not complete) else "  <-- MISMATCH"
                    if complete and e != g:
                        mismatches += 1
                    print(f"   ΔG {dG:+.2f}   collected {g:>10,}   genome {e:>10,}{flag}")
                if complete:
                    print(f"   {'✅ all ΔG counts match' if mismatches == 0 else f'❌ {mismatches} mismatched'}")
            except Exception as ex:
                print(f"   skipped ({ex})")
                complete, mismatches = False, 0

            print("\n🔍 Data quality:")
            checks = {}
            for col in ("value", "record_id", "step_index", "pair"):
                cur.execute(f"SELECT COUNT(*) FROM {TABLE} WHERE {col} IS NULL;")
                checks[f"NULL {col}"] = cur.fetchone()[0]
            cur.execute(f"SELECT COUNT(*) FROM {TABLE} WHERE value > 0;")
            checks["positive ΔG (should be 0)"] = cur.fetchone()[0]
            cur.execute(f"SELECT count_distinct(step_index) FROM {TABLE};")
            distinct = cur.fetchone()[0]
            checks["duplicate step_index"] = total - distinct
            checks["repeated chain_index"] = n - n_distinct
            checks["chain gaps"] = (hi - lo + 1) - n_distinct
            checks["fork imbalance (>1)"] = max(0, fork_imbalance - 1)

            for k, v in checks.items():
                print(f"   {k}: {v:,}")

            if all(v == 0 for v in checks.values()):
                print("✅ Data quality PASSED")
            else:
                print("⚠️  Data quality issues detected")

        conn.close()

    except Exception as e:
        print(f"❌ Data validation failed: {e}")


if __name__ == "__main__":
    validate_genome_data()
