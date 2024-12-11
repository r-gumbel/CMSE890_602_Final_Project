#Snakefile
import os

# Default rule
rule all:
    input:
        "test_tdhf.slurm"

'''configfile' contains dictionaries 
to populate both our .slurm script as 
well as our tdhf3d.inp'''

configfile: 'config.yaml'
slurm_config=config['sconfig']

rule generate_slurm_script:
    output:
        "test_tdhf.slurm"
    run:
        # Generate module load commands
        module_cmds = '\n'.join([f'module load {module}' for module in config.get('modules', [])])
        with open(output[0], 'w') as f:
            f.write(f'''#!/bin/bash
#SBATCH --time={slurm_config['time_limit']}  # limit of wall clock time          
#SBATCH --ntasks={slurm_config['ntasks']}    # number of tasks, i.e. nodes that you require (same as -n)      
#SBATCH --nodes={slurm_config['nodes']}      
#SBATCH --mem-per-cpu={slurm_config['mem_per_cpu']} # memory required per allocated CPU (core) in bytes
#SBATCH --cpus-per-task={slurm_config['cpus-per-task']}          
#SBATCH --job-name={slurm_config['job-name']}
#SBATCH -A {slurm_config['A']}''')
            f.write('\n')
            f.write('''########## Command Lines to Run ##########

export OMP_STACKSIZE=512MB
export OMP_NUM_THREADS=80
ulimit -s unlimited
module purge\n''')
            f.write('\n')
            f.write(f'''# Load required modules
{module_cmds}''')
            f.write('\n\n')
            f.write('''# cp -r /mnt/research/ascsn/Richard/VU-TDHF3D.template /mnt/scratch/gumbelri/static_runs/inverseQF/Ag107_S
 cd /mnt/scratch/gumbelri/static_runs/inverseQF/Ag107_S
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
  107.0D0 12.0D0                                fmass(if) if=1,nof
  47.0D0 6.0D0                                 fcharg(if) if=1,nof
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
  10000      0.400D0                                nt,dt" > /mnt/scratch/gumbelri/static_runs/inverseQF/Ag107_S/run/tdhf3d.inp

 cd /mnt/scratch/gumbelri/static_runs/inverseQF/Ag107_S/run
# srun run 
# scontrol show job $SLURM_JOB_ID     ### write job information to SLURM output file
# js -j $SLURM_JOB_ID                 ### write resource usage to SLURM output file (powetools command)
# ''')
