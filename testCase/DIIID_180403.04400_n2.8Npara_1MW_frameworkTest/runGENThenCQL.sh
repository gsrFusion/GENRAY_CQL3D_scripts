  # expects the target directory to be passed in in which to find the .pbs file 
  genray_id=$(sbatch --parsable $1/genr_sam.pbs)
  cql_id=$(sbatch --parsable --dependency=afterok:$genray_id $1/cql.pbs)