#Snakefile
"""
This is our workflow that generates the slurm scripts for individual
structure calculations of nuclei in their ground state.
It can be used in tandem with the main Snakefile, which runs a larger
workflow of job creation, submission, monitoring, and benchmarking.
Or this .smk can be ran as a standalone slurm script generator for 
a desired nucleus by entering into the command line either:

`snakemake -f -s tdhf.smk`

To use defaults in the config.yaml file.  Or:

`snakemmake -f -s tdhf.smk --config nucleus.A=[Atomic Mass] nucleus.Z=[Atomic Number]`

To override the neutron/proton configuration on file. 
"""
import os
import sys

# # Check if this is being run directly or included
is_run_directly = sys.argv[0].endswith(".smk")

# # Load config with fallback values
configfile: "config.yaml"

# # Set defaults if not provided by main workflow
# if "nucleus" not in config:
#     config["nucleus"] = {"A": 40, "Z": 20}
# if "skyrme" not in config:
#     config["skyrme"] = "SLy4dL"
# if "sconfig" not in config:
#     config["sconfig"] = {
#         "time_limit": "02:00:00",
#         "ntasks": 1,
#         "nodes": 1, 
#         "mem_per_cpu": "2G",
#         "cpus-per-task": 8,
#         "A": "ptg"
#     }


# Define global variables
working_dir = os.getcwd()
run_dir = f'A_{config["nucleus"]["A"]}_Z_{config["nucleus"]["Z"]}_{config["skyrme"]}'
run_id = f'A_{config["nucleus"]["A"]}_Z_{config["nucleus"]["Z"]}_{config["skyrme"]}'


# Create SLURM log directory structure
slurm_log_dir = f"{working_dir}/{run_dir}/slurm_logs"
# os.makedirs(slurm_log_dir, exist_ok=True)

# Create dynamic job name
job_name = f"TDHF_A__{config['nucleus']['A']}_Z_{config['nucleus']['Z']}"

# # A conditionally defined 'rule all' for when tdhf.smk is ran directly
# if is_run_directly:
#     rule all:   
#         input:
#             f"test_tdhf_{run_id}.slurm"

# if is_run_directly==False:
#     os.makedirs(slurm_log_dir, exist_ok=True)

slurm_config=config['sconfig']

rule generate_slurm_script:
    output:
        f"test_tdhf_{run_id}.slurm"
    run:
        # Generate module load commands
        module_cmds = '\n'.join([f'module load {module}' for module in config.get('modules', [])])
        with open(output[0], 'w') as f:
            f.write(f'''#!/bin/bash --login
#SBATCH --time={slurm_config['time_limit']}  # limit of wall clock time          
#SBATCH --ntasks={slurm_config['ntasks']}    # number of tasks, i.e. nodes that you require (same as -n)      
#SBATCH --nodes={slurm_config['nodes']}      
#SBATCH --mem-per-cpu={slurm_config['mem_per_cpu']} # memory required per allocated CPU (core) in bytes
#SBATCH --cpus-per-task={slurm_config['cpus-per-task']}          
#SBATCH --job-name={job_name}
#SBATCH --output={slurm_log_dir}/slurm-%j.out
#SBATCH --error={slurm_log_dir}/slurm-%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user={slurm_config.get('mail_user', '')}
#!SBATCH -A {slurm_config['A']}''')
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
            f.write(f"""cp -r {working_dir}/VU-TDHF3D.template {working_dir}/{run_dir}
cd {working_dir}/{run_dir}
            """)
            f.write('\n')
            f.write('''pwd
./clean
./build.ifc_omp_hpcc''')

            f.write('''
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
00.0   00.0                                   euler_gamma(if)\n''')
            f.write(f'''{config['nucleus']['A']}.0D0 12.0D0                                fmass(if) if=1,nof
{config['nucleus']['Z']}.0D0 6.0D0                                 fcharg(if) if=1,nof\n''')
            f.write('''3.0 3.0                                       radinf(1,if) if=1,nof
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
20 20                                         nextra_n(if), if=1,nof
20 20                                         nextra_p(if), if=1,nof
1                                             ifixb
15.0  20.0   0.0                             ecm,rsep,xb
12     1.0D-17                                mxp,terr
10000      0.400D0                                nt,dt" > ''') 
            f.write(f'''{working_dir}/{run_dir}/run/tdhf3d.inp \n

cd {working_dir}/{run_dir}/run/''')
            f.write('''\npwd\n
srun run 
scontrol show job $SLURM_JOB_ID     ### write job information to SLURM output file
js -j $SLURM_JOB_ID                 ### write resource usage to SLURM output file (powertools command)''')