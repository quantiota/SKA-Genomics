

WITH base_data AS (
  SELECT
    r.step_index,
    r.chain_index,
    r.entropy,
    s.pair,
    s.transition_index,
    LAG(r.entropy) OVER (ORDER BY r.step_index) AS prev_entropy
  FROM genome_ska_results r
  JOIN genome_steps s ON (record_id, fork, step_index)
  WHERE r.record_id = 'NC_000913.3'
    AND r.fork      = 'fork1'
    AND r.entropy IS NOT NULL
),
with_probability AS (
  SELECT
    step_index,
    chain_index,
    entropy,
    pair,
    transition_index,
    CASE
      WHEN entropy != 0 AND prev_entropy IS NOT NULL
      THEN EXP(-ABS((entropy - prev_entropy) / entropy))
      ELSE NULL
    END AS P
  FROM base_data
  WHERE prev_entropy IS NOT NULL
)
SELECT
  step_index,
  chain_index,
  P,
  transition_index,
  concat(left(pair, 1), '->', right(pair, 1)) AS transition_name
FROM with_probability
ORDER BY step_index;


