

SELECT
  step_index,
  entropy,
  fork
FROM genome_ska_results
WHERE record_id = 'NC_000913.3'
  AND entropy IS NOT NULL
ORDER BY step_index;


