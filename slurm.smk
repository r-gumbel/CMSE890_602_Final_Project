#Snakefile
import os
# Default rule

rule all:
    input:
        "run_routine.slurm"

# Slurm configuration dictionary
slurm_config = {
    'job_name': 'fortran_hello',
    'output_log': 'hello_output_%j.txt',
    'error_log': 'hello_error_%j.txt',
    'time_limit': '00:05:00', 
    'ntasks': 1,
    'mem': '1G',
    'partition': 'short',
    'modules': [
        # List of modules to load, if any
        # 'gcc/10.2.0'
    ]
}

# Configuration
EXECUTABLE_NAME = 'hello_world'

# Rule to generate some test fortran source code
rule create_fortran_source:
    output:
        "src/main.f90"
    run:
        # Ensure the src directory exists
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)

        # Write the fortran source code 
        with open(output[0], 'w') as f:
            f.write('''program hello_world
    print *, "Hello World!"
end program hello_world''')


#Rule to compile the fortran code
rule compile_fortran:   
    input:
        "src/main.f90"
    output:
        EXECUTABLE_NAME
    shell:
        "gfortran {input} -o {output}"


rule generate_slurm_script:
    input:
        exe=EXECUTABLE_NAME
    output:
        "run_routine.slurm"
    run:
        # Generate module load commands
        module_cmds = '\n'.join([f'module load {module}' for module in slurm_config.get('modules', [])])

        with open(output[0], 'w') as f:
            f.write(f'''#!/bin/bash
#SBATCH --job-name={slurm_config['job_name']}
#SBATCH --output={slurm_config['output_log']}
#SBATCH --error={slurm_config['error_log']}
#SBATCH --time={slurm_config['time_limit']}
#SBATCH --ntasks={slurm_config['ntasks']}
#SBATCH --mem={slurm_config['mem']}
#SBATCH --partition={slurm_config['partition']}

# Print some job information
echo "Job started on $(data)"
echo "Running on host $(hostname)"
echo "Job ID: $SLURM_JOB_ID"

# Load required modules
{module_cmds}

# Run the executable
./{input.exe}

# Print the completion time
echo "Job completed on $(date)"
''')
