#!/bin/bash --login
#SBATCH --time=0:10:00  # limit of wall clock time          
#SBATCH --ntasks=1    # number of tasks, i.e. nodes that you require (same as -n)      
#SBATCH --nodes=1      
#SBATCH --mem-per-cpu=8G # memory required per allocated CPU (core) in bytes
#SBATCH --cpus-per-task=12          
#SBATCH --job-name=TDHF_test
#SBATCH -A ptg
########## Command Lines to Run ##########

export OMP_STACKSIZE=512MB
export OMP_NUM_THREADS=80
ulimit -s unlimited
module purge

# Load required modules
module load powertools
module load intel/2023a
module load FFTW

cp -r /mnt/home/gumbel/CMSE890_602_Final_Project/VU-TDHF3D.template /mnt/home/gumbelri/CMSE890_602_Final_Project/test
cd /mnt/home/gumbel/CSME890_602_Final*/test
pwd
./clean
./build.ifc_omp_hpcc

echo "SLy4dL
  1                                             nof
  0                                             imode
   0.0    0.0                                   centf(1,if) if=1,nof
   0.0    0.0                                   centf(2,if) if=1,nof
   0.0    0.0                                   centf(3,if) if=1,nof
   0.0    0.0                                   boostf(1,if) if=1,nof
   0.0    0.0                                   boostf(2,if) if=1,nof
   0.0    0.0                                   boostf(3,if) if=1,nof
  00.0   00.0                                   euler_alpha(if)
  00.0   00.0                                   euler_beta(if)
  00.0   00.0                                   euler_gamma(if)
  4.0D0 12.0D0                                fmass(if) if=1,nof
  2.0D0 6.0D0                                 fcharg(if) if=1,nof
  3.0 3.0                                       radinf(1,if) if=1,nof
  3.0 3.0                                       radinf(2,if) if=1,nof
  3.0 3.0                                       radinf(3,if) if=1,nof
  0.02   35.0                                  x0dmp,e0dmp
   7500  2500  6.0D-4 7.0D-2                   itrbx,mtrbx,serr,derr
   1   20   00  20                               iprint,mprint,mplots,mplott
   0 0 200 200                              irest,irwgs,mrests,mrestt
   0 1                                          ifixcm,icoul
   0 0.000 20  1.90 5e-5                               iconstr,q2in,mconstr,c0=(1.9),d0=(5e-5)
   0 0 1 1                                      itimrev, itimrevs, iodds, hdiag
   0 0                                          ipairf(1,if) if=1,nof
   0 0                                          ipairf(2,if) if=1,nof
   0 2000                                       nexadd, nexiter
   26 20                                         nextra_n(if), if=1,nof
   16 20                                         nextra_p(if), if=1,nof
  1                                             ifixb
  15.0  20.0   0.0                             ecm,rsep,xb
  12     1.0D-17                                mxp,terr
  10000      0.400D0                                nt,dt" > /mnt/home/gumbelri/CMSE890_602_Final*/test/run/tdhf3d.inp

cd /mnt/home/gumbelri/CMSE890_602_Final*/test/run
pwd
# srun run 
# scontrol show job $SLURM_JOB_ID     ### write job information to SLURM output file
# js -j $SLURM_JOB_ID                 ### write resource usage to SLURM output file (powetools command)
# 