process PLOT_GTCHECK_SCOREDIST
{
  input:
    path(pool_scores_csv)

  script:
  """
    Rscript --vanilla ${projectDir}/../bin/gtcheck_scoredist.R ${pool_scores_csv}
  """
}
