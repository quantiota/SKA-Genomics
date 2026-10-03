"""Configuration for genomic data validation."""
import os

# QuestDB ILP (line protocol) — ingestion path, 9009
ILP_CONFIG = {
    'host': os.getenv('QDB_ILP_HOST', os.getenv('QDB_PG_HOST', 'localhost')),
    'port': int(os.getenv('QDB_ILP_PORT', '9009')),
}

# QuestDB pg-wire — queries and DDL, 8812
QDB_CONFIG = {
    'dbname': os.getenv('QDB_PG_NAME', 'qdb'),
    'user': os.getenv('QDB_PG_USER', 'admin'),
    'password': os.getenv('QDB_PG_PASSWORD', 'quest'),
    'host': os.getenv('QDB_PG_HOST', 'localhost'),
    'port': os.getenv('QDB_PG_PORT', '8812'),
}

# Genomic Stream Configuration
GENOME_CONFIG = {
    'fasta': os.getenv('GENOME_FASTA',
                       'data/GCF_000005845.2_ASM584v2_genomic.fna'),
    'record_id': 'NC_000913.3',      # the FASTA record identifier
    'organism': 'E. coli K-12 MG1655',
    # A chain has no native sampling rate, so one is IMPOSED to make the stream
    # real-time. Timestamps are then real wall-clock, and delta_t downstream is
    # MEASURED from them as (t_k - t_{k-1}).total_seconds() rather than assumed.
    # 1,000 steps/s = the E. coli replication fork rate (DNA Pol III), so
    # delta_t = 1 ms is the real interval at which the cell builds the chain.
    # Cross-check: two forks x 2.32 Mb / 2,400 s (the ~40 min C-period) ~ 967 bp/s.
    # One single-pass walk of 4,641,651 steps therefore takes ~77 min.
    'rate_hz': 1000.0,               # steps/s — biological replication rate
    'packet_steps': 1000,            # one packet per second
    'max_steps': None,               # None = whole chain
    # Emission order. 'linear' walks position 1 -> end, which is NOT the order
    # the cell reads the chain. 'replication' sends two forks outward from oriC.
    'order': 'linear',
    'origin': 3925744,               # E. coli oriC, approx. chain index
}

# Logging Configuration
LOG_CONFIG = {
    'level': 'INFO',
    'format': '%(asctime)s %(levelname)s %(message)s',
    'file': 'logs/genome_validation.log',
}
