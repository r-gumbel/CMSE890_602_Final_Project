#Snakefile
"""
This is our workflow that generates the slurm scripts for individual
structure calculations of nuclei in their ground state.
"""
import os
from pathlib import Path

# Define variable for working directory (one level up)
working_dir = str(Path(os.getcwd()).parent)
# Define run directory, following naming convetion of main workflow
run_dir = f'A_{config["nucleus"]["A"]}_Z_{config["nucleus"]["Z"]}_{config["skyrme"]}'

# Create dynamic job name (fixed indentation)
job_name = f"TDHF_A_{config['nucleus']['A']}_Z_{config['nucleus']['Z']}"
 
# Default rule
rule script_maker:
    input:
        f"../tdhf_{run_id}.slurm"     # run_id comes from main workflow

slurm_config=config['sconfig']

rule generate_slurm_script:
    output:
        f"../tdhf_{run_id}.slurm"
    run:
        # Generate module load commands
        module_cmds = '\n'.join([f'module load {module}' for module in config.get('modules', [])])
        
        # Get email address from config
        email_address = slurm_config.get('mail_user', '')
        
        with open(output[0], 'w') as f:
            f.write(f'''#!/bin/bash --login
#SBATCH --time={slurm_config['time_limit']}  # limit of wall clock time          
#SBATCH --ntasks={slurm_config['ntasks']}    # number of tasks, i.e. nodes that you require (same as -n)      
#SBATCH --nodes={slurm_config['nodes']}      
#SBATCH --mem-per-cpu={slurm_config['mem_per_cpu']} # memory required per allocated CPU (core) in bytes
#SBATCH --cpus-per-task={slurm_config['cpus-per-task']}          
#SBATCH --job-name={job_name}
#SBATCH --mail-type=END
#SBATCH --mail-user=''')
            # Write email address directly
            f.write(email_address)
            # Continue with rest of SLURM directives
            f.write(f'''
#!SBATCH -A {slurm_config['A']}
            
########## Command Lines to Run ##########

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
2500  2500  6.0D-3 7.0D-2                   itrbx,mtrbx,serr,derr
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
js -j $SLURM_JOB_ID                 ### write resource usage to SLURM output file (powertools command)

# Wait for the report to be generated
sleep 60  # Give some time for report generation

# Get absolute paths
REPORT_PATH="${PWD%/*}/logs/${SLURM_JOB_NAME}/report.txt"
echo "Looking for report at: $REPORT_PATH"

# Check if report exists and email it
if [ -f "$REPORT_PATH" ]; then
    if [ ! -z "$SLURM_JOB_MAIL_USER" ]; then
        echo "Sending report to $SLURM_JOB_MAIL_USER"
        mail -s "TDHF Job $SLURM_JOB_ID Report" "$SLURM_JOB_MAIL_USER" < "$REPORT_PATH"
    else
        echo "No email address configured - report will not be sent"
    fi
else
    echo "Report file not found at $REPORT_PATH"
    if [ ! -z "$SLURM_JOB_MAIL_USER" ]; then
        echo "Report file not found at $REPORT_PATH" | mail -s "TDHF Job $SLURM_JOB_ID - Report Missing" "$SLURM_JOB_MAIL_USER"
    else
        echo "No email address configured - error report will not be sent"
    fi
fi
''')