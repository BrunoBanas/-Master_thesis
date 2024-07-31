#!/bin/bash
#SBATCH -J vasp_job_1100
#SBATCH -p c8
#SBATCH -N 2
#SBATCH --ntasks-per-node=40
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err


# Load modules
module purge
module load intel/2021.1.2.266 openmpi4/4.0.4 vasp/6.2.0

# Function to check XDATCAR and rerun job if necessary
check_and_rerun() {
    if grep -q "Direct configuration=     2" XDATCAR; then
        echo "Line found in XDATCAR. Copying CONTCAR to POSCAR and rerunning the job."
        cp CONTCAR POSCAR
        # Rerun the job
        mpirun -n 40 vasp_std
    else
        echo "Line not found in XDATCAR. Job completed."
    fi
}

# Run job
mpirun -n 80 vasp_std

# Check XDATCAR and rerun if necessary
# check_and_rerun
