#!/bin/bash
#SBATCH --job-name=fortran_hello
#SBATCH --output=hello_output_%j.txt
#SBATCH --error=hello_error_%j.txt
#SBATCH --time=00:05:00
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --partition=short

./hello_world
