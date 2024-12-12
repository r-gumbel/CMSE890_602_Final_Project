# Main snakefile 
import os
import subprocess
from datetime import datetime
from pathlib import Path 

# Include the TDHF workflow
include: "tdhf.smk"

# Create basic log structure
Path("logs").mkdir(exist_ok=True)

def submit_and_get_jobid(script_path):
    """Submit SLURM job and return job ID"""
    result = subprocess.run(['sbatch', script_path],
                            capture_output=True, text=True)
    job_id = result.stdout.strip().split()[-1]
    return job_id

rule all:
    input:
        "logs/tdhf_status.txt"

rule run_tdhf:
    input:
        slurm_script = "test_tdhf.slurm"
    output:
        status = "logs/tdhf_status.txt"
    run:
        # Submit the job
        job_id = submit_and_get_jobid(input.slurm_script)

        # Log the submission
        with open(output.status, 'w') as f:
            f.write(f"""TDHF Job Submitted
Timestamp: {datetime.now()}
Job ID: {job_id}
Script: {input.slurm_script}
""")

        print(f"Submitted TDHF job {job_id}")
        print(f"Check status with: squeue -j {job_id}")
        print(f"Check logs at: logs/tdhf_status.txt")