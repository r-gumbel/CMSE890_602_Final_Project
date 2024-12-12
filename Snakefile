# Main snakefile 
import os
import subprocess
from datetime import datetime
from pathlib import Path 

"""
This is the main workflow, which coordinates all the
processes involved, including separate workflows:
    (i.) tdhf.smk (creates a slurm script to run VU-TDHF)
    (ii.) report.smk (looks over static HF results and prepares a report)

This main Snakefile workflow will be executed 
by typing:
    `snakemake -j1 --config A=48 Z=20`
"""

# configuration file
configfile: 'config.yaml'

# Debug print initial config
print("Initial config state:")
print(config)

# Handle command line overrides for A and Z
if 'A' in config:
    config['nucleus']['A'] = int(config.pop('A'))
if 'Z' in config:
    config['nucleus']['Z'] = int(config.pop('Z'))

# Set defaults if not provided
if 'nucleus' not in config:
    config['nucleus'] = {}
    
config['nucleus'].setdefault('A', 40)
config['nucleus'].setdefault('Z', 20)
config['skyrme'] = config.get('skyrme', 'SLy4dL')

# Debug print final configuration
print("\nFinal configuration:")
print(f"A = {config['nucleus']['A']}")
print(f"Z = {config['nucleus']['Z']}")
print(f"Skyrme = {config['skyrme']}")

# Generate run ID automatically from parameters
run_id = f"A_{config['nucleus']['A']}_Z_{config['nucleus']['Z']}_{config['skyrme']}"

# Include the TDHF workflow
include: "tdhf.smk"

# Create basic log structure
Path("logs/{run_id}").mkdir(parents=True, exist_ok=True)

def submit_and_get_jobid(script_path):
    """Submit SLURM job and return job ID"""
    result = subprocess.run(['sbatch', script_path],
                          capture_output=True, text=True)
    job_id = result.stdout.strip().split()[-1]
    return job_id

rule all:
    input:
        f"logs/{run_id}/tdhf_status.txt"

rule run_tdhf:
    input:
        slurm_script = f"test_tdhf_{run_id}.slurm"
    output:
        status = f"logs/{run_id}/tdhf_status.txt"
    run:
        # Submit the job
        job_id = submit_and_get_jobid(input.slurm_script)
        
        # Log the submission
        with open(output.status, 'w') as f:
            f.write(f"""TDHF Job Submitted
Timestamp: {datetime.now()}
Job ID: {job_id}
Run ID: {run_id}
Script: {input.slurm_script}
Nucleus: A={config['nucleus']['A']}, Z={config['nucleus']['Z']}
Skyrme: {config['skyrme']}
""")

        print(f"Submitted TDHF job {job_id}")
        print(f"Check status with: squeue -j {job_id}")
        print(f"Check logs at: {output.status}")