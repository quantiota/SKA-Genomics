-- Genomic transition probability bands
-- Datasource: QuestDB plugin  ·  Tables: genome_ska_results (SKA output)
--                                        genome_steps       (source chain)
--
-- Panel: XY Chart
--   x = step_index          (emission order; fork1 holds the odd indices)
--   y = P                   transition probability in (0, 1]
--   Draw style: POINTS.
--   Color: transform "Partition by values" on `transition_name`
--          -> 16 colored, legended bands.
--
-- The market template derives its regime from the sign of the price change,
-- giving 3 states and 3x3 = 9 transitions. Here the state is the base itself,
-- so there are 4 states and 4x4 = 16 transitions -- and the transition is not
-- derived at all: genome_steps already stores it as `pair` ('AT') with its
-- `transition_index` (0-15, letter order A, T, G, C). The LAG(regime) stage of
-- the template is therefore unnecessary; only LAG(entropy) is kept.
--
--   P = exp(-|(H_k - H_{k-1}) / H_k|)
--
-- P -> 1 when entropy barely moved across the step (the transition was
-- expected), P -> 0 on a large relative jump (the transition was surprising).

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


-- Band summary -- mean P per transition over the learned window.
-- Panel: Bar chart, x = transition_name, y = mean_P.
--
-- WITH base_data AS (
--   SELECT r.step_index, r.entropy, s.pair,
--          LAG(r.entropy) OVER (ORDER BY r.step_index) AS prev_entropy
--   FROM genome_ska_results r
--   JOIN genome_steps s ON (record_id, fork, step_index)
--   WHERE r.record_id = 'NC_000913.3' AND r.fork = 'fork1'
--     AND r.entropy IS NOT NULL
-- )
-- SELECT
--   concat(left(pair, 1), '->', right(pair, 1)) AS transition_name,
--   count()                                     AS steps,
--   avg(EXP(-ABS((entropy - prev_entropy) / entropy)))  AS mean_P,
--   min(EXP(-ABS((entropy - prev_entropy) / entropy)))  AS min_P,
--   max(EXP(-ABS((entropy - prev_entropy) / entropy)))  AS max_P
-- FROM base_data
-- WHERE prev_entropy IS NOT NULL
-- GROUP BY pair
-- ORDER BY mean_P DESC;
