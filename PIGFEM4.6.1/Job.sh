#!/bin/bash
#SBATCH -p normal_q
#SBATCH --ntasks=1
#SBATCH -A shakiba_cfrp
#SBATCH -t 200:00:00
#SBATCH --mem-per-cpu=50G

DIR_ADRS=$(pwd)
export LD_LIBRARY_PATH=${DIR_ADRS}/petsc-3.7.3/arch_external_opt/lib:${DIR_ADRS}/gsl-2.1/install/lib
module load gcc/6.1.0

ulimit -s unlimited

./Main -in "Input/3_531_fibers.csv" -sensitivity true