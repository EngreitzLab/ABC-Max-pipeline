#!/bin/bash
#SBATCH --job-name=ABC-Max             # Job name
#SBATCH --partition=akundaje                      # Partition name
#SBATCH --time=2-00:00                    # Runtime in D-HH:MM format
#SBATCH --nodes=1                         # Number of nodes (keep at 1)
#SBATCH --ntasks=1                        # Number of tasks per node (keep at 1)
#SBATCH --cpus-per-task=1                 # CPU cores requested per task (change for threaded jobs)
#SBATCH --mem=12G                        # Memory needed per node (total)
#SBATCH --error=jobid_%j.err              # File to which STDERR will be written, including job ID
#SBATCH --output=jobid_%j.out             # File to which STDOUT will be written, including job ID

#source activate abc-max

snakemake --snakefile ABC-Max.snakefile --configfile ../config/ABC-Max.example.json -j 1 --keep-target-files --rerun-incomplete --cluster "sbatch -n 1 -c 12 --mem 57G -t 2-00:00"
