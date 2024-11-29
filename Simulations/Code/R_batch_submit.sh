#!/bin/bash
#SBATCH -J WRPC_scen21_        # Job name for the array
#SBATCH -o WRPC_scen21_%A.out  # Shared standard output with job ID
#SBATCH -e WRPC_scen21_%A.err  # Shared standard error with job ID
#SBATCH -p shared      # Partition to submit to
#SBATCH -c 1	       # Number of cores (for parallelization)
#SBATCH -N 1           # Number of nodes
#SBATCH -t 0-4:00:00  # Runtime (D-HH:MM:SS)
#SBATCH --mem=9000     # Memory request
#SBATCH --mail-type=BEGIN,END,FAIL  # Mail notifications
#SBATCH --mail-user=stephaniewu@fas.harvard.edu  # Account to email

module load R/4.2.2-fasrc01 gcc/10.2.0-fasrc01
export R_LIBS_USER=$HOME/apps/R_4.2.2:$R_LIBS_USER
scenario=21
Rscript /n/netscratch/stephenson_lab/Lab/stephwu18/WRPC/Simulations/Code/run_sims.R ${scenario} ${SLURM_ARRAY_TASK_ID}
