  # expects the target directory to be passed in in which to find the .pbs file 
  cql_id=$(sbatch --parsable $1/cql.pbs)