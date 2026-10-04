"""Simulated genomic stream into QuestDB — raw parameters only.

A genomic chain has no native sampling rate, so one is IMPOSED. Steps are
emitted in packets paced against the wall clock, and each row is stamped on a
FIXED GRID anchored at the run's start time:

    timestamp(k) = t_start + k / rate_hz

These are ASSIGNED timestamps, not observations. delta_t downstream will
therefore return exactly 1/rate_hz with no jitter — it recovers the configured
rate rather than measuring anything. That determinism is deliberate and is the
right property for analysis, but it must not be mistaken for a measurement. To
obtain real jitter the stream would have to stamp each row with time.time_ns()
at the moment of emission, one row at a time.

Emission order is selectable:

    linear       position 1 -> end. A plain walk. NOT how the cell reads it.
    replication  two forks leaving oriC in opposite directions around the
                 circular chromosome, interleaved, meeting at the terminus.

No returns, no entropy, no learning. Raw values only.
"""
import argparse
import datetime
import logging
import os
import signal
import time
import traceback

import socket

import psycopg2
import psycopg2.pool

from config import QDB_CONFIG, ILP_CONFIG, GENOME_CONFIG, LOG_CONFIG

from genome_chain import GenomeChain

os.makedirs('logs', exist_ok=True)
logging.basicConfig(level=getattr(logging, LOG_CONFIG['level']),
                    format=LOG_CONFIG['format'],
                    filename=LOG_CONFIG['file'], filemode='a')
console = logging.StreamHandler()
console.setLevel(logging.INFO)
console.setFormatter(logging.Formatter(LOG_CONFIG['format']))
logging.getLogger('').addHandler(console)

connection_pool = psycopg2.pool.SimpleConnectionPool(1, 10, **QDB_CONFIG)

TABLE = "genome_steps"
FORK_TAG = {0: "linear", 1: "fork1", 2: "fork2"}

CREATE_SQL = f"""
CREATE TABLE {TABLE} (
    record_id SYMBOL,
    organism SYMBOL,
    step_index LONG,
    fork SYMBOL,
    chain_index LONG,
    position LONG,
    pair SYMBOL,
    transition_index LONG,
    value DOUBLE,
    twist DOUBLE,
    tilt DOUBLE,
    level DOUBLE,
    rate_hz DOUBLE,
    chain_length LONG,
    total_steps LONG,
    timestamp TIMESTAMP
) TIMESTAMP(timestamp) PARTITION BY DAY;
"""


def ilp_escape(v: str) -> str:
    """Escape an ILP tag value: spaces, commas and equals must be backslashed."""
    return str(v).replace("\\", "\\\\").replace(" ", "\\ ").replace(",", "\\,").replace("=", "\\=")


def ilp_connect():
    """Open the ILP ingestion socket (port 9009). ~600k rows/s vs ~3.5k on pg-wire."""
    return socket.create_connection((ILP_CONFIG["host"], ILP_CONFIG["port"]), timeout=30)


def create_questdb_table(recreate=False):
    """Create the table. Refuses to drop a non-empty one unless recreate=True:
    an unconditional DROP here will silently destroy a collection that another
    process is still writing."""
    conn = connection_pool.getconn()
    try:
        with conn.cursor() as cur:
            existing = 0
            try:
                cur.execute(f"SELECT count() FROM {TABLE};")
                existing = cur.fetchone()[0]
            except Exception:
                conn.rollback()                 # table does not exist yet
            if existing and not recreate:
                raise SystemExit(
                    f"{TABLE} already holds {existing:,} rows. Pass --recreate to "
                    f"drop it, or use a different table. Refusing to destroy data.")
            if existing:
                logging.warning(f"--recreate: dropping {existing:,} existing rows")
            cur.execute(f"DROP TABLE IF EXISTS {TABLE};")
            cur.execute(CREATE_SQL)
            conn.commit()
            logging.info(f"Clean {TABLE} table created")
    except Exception as e:
        logging.error(f"Error creating QuestDB table: {e}")
        logging.error(traceback.format_exc())
    finally:
        connection_pool.putconn(conn)


def send_packet(sock, lines):
    """Ship one packet over ILP. Raises on failure — a dropped packet is a gap
    in the chain, and TCP ILP does not report rejected rows, so continuing
    would silently corrupt the collection."""
    try:
        sock.sendall("".join(lines).encode())
    except Exception as e:
        logging.error(f"Error sending genomic packet over ILP: {e}")
        logging.error(traceback.format_exc())
        raise


def run_chain_stream(fasta, record_id, organism, rate_hz, packet_steps,
                     max_steps, order="linear", origin=None,
                     leading_strand=True):
    """Emit the chain in packets, paced against the wall clock, over ILP."""
    chain = GenomeChain(fasta, max_steps=max_steps, order=order, origin=origin,
                        leading_strand=leading_strand)
    chain_length, total_steps = chain.chain_length, chain.total_steps

    # (2) the record ID comes from the FASTA header; a mismatch is fatal
    accession = chain.accession
    if record_id and record_id != accession:
        raise SystemExit(f"record_id mismatch: config says {record_id!r}, "
                         f"FASTA header says {accession!r}. Fix GENOME_CONFIG "
                         f"['record_id'] or point at the right file.")
    record_id = accession

    # (4) duration reflects what will actually be emitted
    n_emit = total_steps if max_steps is None else min(max_steps, total_steps)
    dt = 1.0 / rate_hz
    logging.info(f"Streaming {record_id} | {chain_length:,} bp | {total_steps:,} steps | "
                 f"{rate_hz:,.0f} Hz | packets of {packet_steps:,} | order={order}")
    n_forks = 2 if order == "replication" else 1
    logging.info(f"Emitting {n_emit:,} steps -> projected duration "
                 f"{n_emit / rate_hz / 60:.2f} min")
    if n_forks == 2:
        logging.info(f"{rate_hz:,.0f} Hz TOTAL is shared between 2 forks = "
                     f"{rate_hz / 2:,.0f} bp/s per fork. Use --rate 2000 for the "
                     f"~1,000 bp/s biological fork rate (~39 min, matching the "
                     f"~40 min C-period).")

    rec, org = ilp_escape(record_id), ilp_escape(organism)
    sock = ilp_connect()
    t_start = time.time()
    t0_ns = int(t_start * 1e9)

    lines, count, packets = [], 0, 0
    try:
        for e in chain.events():
            # (1) ASSIGNED timestamp on a fixed grid. (6c) computed from the
            # exact fraction, so a non-integer rate does not accumulate drift.
            ts_ns = t0_ns + round(count * 1e9 / rate_hz)
            lines.append(
                f"{TABLE},record_id={rec},organism={org},fork={FORK_TAG[e['fork']]},"
                f"pair={e['pair']} "
                f"step_index={e['step']}i,chain_index={e['chain_index']}i,"
                f"position={e['position']}i,transition_index={int(e['index'])}i,"
                f"value={e['value']},"
                f"twist={e['twist']},tilt={e['tilt']},level={e['level']},"
                f"rate_hz={rate_hz},chain_length={chain_length}i,"
                f"total_steps={total_steps}i {ts_ns}\n")
            count += 1

            if len(lines) >= packet_steps:
                lag = (t_start + count * dt) - time.time()
                if lag > 0:
                    time.sleep(lag)
                send_packet(sock, lines)
                lines = []
                packets += 1
                if packets % 60 == 0:
                    el = time.time() - t_start
                    logging.info(f"Collected {count:,} raw genomic steps  "
                                 f"({count/el:,.0f} steps/s, {el:.0f}s elapsed)")
    finally:
        # (6b) never drop the tail, including on Ctrl-C
        if lines:
            try:
                send_packet(sock, lines)
                logging.info(f"Flushed final partial packet ({len(lines)} steps)")
            except Exception:
                logging.error("Final packet could not be flushed")
        sock.close()

    elapsed = time.time() - t_start
    logging.info(f"Done: {count:,} steps in {elapsed:.1f}s ({count/elapsed:,.0f} steps/s)")
    return count


def shutdown(signum, frame):
    """Raise, so the stream's finally-block flushes the partial packet."""
    logging.info(f"Received shutdown signal: {signum}")
    raise KeyboardInterrupt


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--rate", type=float, default=GENOME_CONFIG["rate_hz"],
                   help="imposed sampling rate in Hz")
    p.add_argument("--packet", type=int, default=GENOME_CONFIG["packet_steps"])
    p.add_argument("--max-steps", type=int, default=GENOME_CONFIG["max_steps"])
    p.add_argument("--recreate", action="store_true",
                   help="drop an existing non-empty table before collecting")
    p.add_argument("--order", default=GENOME_CONFIG.get("order", "linear"),
                   choices=("linear", "replication"))
    p.add_argument("--origin", type=int, default=GENOME_CONFIG.get("origin"),
                   help="oriC chain index (0-based), used when --order replication")
    p.add_argument("--raw-strand", action="store_true",
                   help="label fork 2 with raw top-strand pairs instead of its "
                        "own leading strand (dG is unaffected either way)")
    a = p.parse_args()

    signal.signal(signal.SIGINT, shutdown)
    signal.signal(signal.SIGTERM, shutdown)

    create_questdb_table(recreate=a.recreate)
    logging.info("Starting CLEAN genomic data stream processor...")
    logging.info(f"Source: {GENOME_CONFIG['fasta']}")
    logging.info("✅ Raw data collection only")
    logging.info("✅ No SKA computations")
    logging.info("✅ Assigned timestamps on a fixed grid (not measured)")
    logging.info(f"✅ ILP ingestion on {ILP_CONFIG['host']}:{ILP_CONFIG['port']}")

    try:
        run_chain_stream(GENOME_CONFIG["fasta"], GENOME_CONFIG["record_id"],
                         GENOME_CONFIG["organism"], a.rate, a.packet, a.max_steps,
                         order=a.order, origin=a.origin,
                         leading_strand=not a.raw_strand)
    except KeyboardInterrupt:
        logging.info("Interrupted — partial packet flushed")
    connection_pool.closeall()
