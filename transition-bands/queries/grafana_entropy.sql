-- SKA entropy evolution along the genomic chain
-- Datasource: QuestDB plugin  ·  Table: genome_ska_results
--
-- Panel: XY Chart
--   x = step_index         (emission order within the series)
--   y = entropy            H(k)
--   Draw style: POINTS, not lines.
--     Consecutive steps belong to different input levels, so joining them
--     paints a solid band. As points the trajectory resolves into a small
--     family of smooth strands — one per discrete input level.
--   Color: transform "Partition by values" on `fork`
--     (one colored, legended series per replication fork).

SELECT
  step_index,
  entropy,
  fork
FROM genome_ska_results
WHERE record_id = 'NC_000913.3'
  AND entropy IS NOT NULL
ORDER BY step_index;


-- Variant: entropy against genomic position rather than emission order.
-- In replication order the two forks run outward from oriC, so this spreads
-- the series across the chromosome instead of stacking them.
--
-- SELECT
--   chain_index,
--   entropy,
--   fork
-- FROM genome_ska_results
-- WHERE record_id = 'NC_000913.3'
--   AND entropy IS NOT NULL
-- ORDER BY chain_index;


-- Variant: the saturation diagnostic. D climbing to 1 means D(1-D) -> 0
-- and the gradient is dead; read alongside the entropy panel.
--
-- SELECT
--   step_index,
--   decision,
--   knowledge,
--   fork
-- FROM genome_ska_results
-- WHERE record_id = 'NC_000913.3'
--   AND entropy IS NOT NULL
-- ORDER BY step_index;
