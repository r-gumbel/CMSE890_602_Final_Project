# Nuclear Structure Calculation Workflow Pseudocode

## Main Workflow (Snakefile)
```plaintext
1. LOAD configuration from config.yaml
2. IF command line arguments provided:
 - OVERRIDE config with command line values for:
   * A (Atomic Mass)
   * Z (Atomic Number)
   * email
   * skyrme force
3. SET default values if not provided:
 - Default A: 40
 - Default Z: 20
 - Default Skyrme Force: SLy4dL
4. GENERATE run_id:
 - Format: "A_{A}_Z_{Z}_{Skyrme}"
5. CREATE logs directory:
 - Path: ../logs/{run_id}
6. INCLUDE secondary workflow files:
 - tdhf.smk
 - report.smk
7. DEFINE helper functions:
 a. check_job_status:
    - USE squeue to verify job completion
 b. wait_for_job:
    - POLL job status every 60 seconds
 c. submit_and_get_jobid:
    - SUBMIT SLURM script
    - RETURN job ID
8. DEFINE primary workflow rules:
 a. all: REQUIRE final report generation
 b. run_tdhf:
    - GENERATE SLURM script
    - SUBMIT job
    - TRACK job status
    - CREATE completion markers
```

## TDHF Script Generation Workflow (tdhf.smk)
```plaintext
1. SET working directory paths:
 - Parent working directory
 - Run-specific directory
2. GENERATE dynamic job name:
 - Based on A, Z values
3. PREPARE SLURM script generation rule:
 a. CONFIGURE SLURM directives:
    - Time limit
    - Resource allocation
    - Job name
    - Email notifications
 b. SET environment variables:
    - OpenMP stack size
    - Thread count
    - System limits
 c. LOAD required modules dynamically
 d. COPY template directory to run location
 e. PREPARE input file (tdhf3d.inp):
    - Nucleus parameters
    - Computational domain settings
    - Convergence parameters
    - Numerical controls
 f. EMBED run commands:
    - Directory navigation
    - Executable calls
    - Job status tracking
    - Report generation
    - Optional email notification
```

## Report Generation Workflow (report.smk)
```plaintext
1. DEFINE data extraction functions:
 a. get_energy:
    - LOCATE most recent output file
    - EXTRACT energy value
    - HANDLE file/parsing errors
 b. get_Q20:
    - EXTRACT quadrupole moment
    - USE grep/awk for parsing
 c. get_iterations:
    - FIND final iteration count
 d. parse_hfb_data:
    - READ reference HFB.data file
    - MATCH nucleus parameters
 e. calculate_relative_error:
    - COMPUTE percentage difference
    - HANDLE edge cases

2. DEFINE report generation rule:
 a. REQUIRE job completion marker
 b. EXTRACT calculation results:
    - Energy
    - Q20
    - Iterations
    - Convergence parameters
 c. PERFORM comparative analysis:
    - Compare with HFB reference data
    - Calculate relative errors
 d. GENERATE comprehensive report:
    - Timestamp
    - Run configuration
    - Calculation results
    - Benchmarking metrics
    - Error analysis
```

## Workflow Execution Sequence
```plaintext
1. Preparation Phase:
 a. LOAD configuration
 b. VALIDATE and PROCESS inputs
 c. GENERATE run identifiers
 d. PREPARE directory structure

2. TDHF Calculation Phase:
 a. GENERATE SLURM script
 b. SUBMIT job to cluster
 c. MONITOR job status
 d. WAIT for completion
 e. CREATE completion markers

3. Analysis Phase:
 a. VERIFY job completion
 b. EXTRACT calculation data
 c. PARSE output files
 d. COMPARE with reference data
 e. GENERATE detailed report

4. Reporting Phase:
 a. CREATE report.txt
 b. OPTIONALLY email results
 c. CLEAN temporary files
```

## Dependency Hierarchy
```plaintext
report.txt
 └── job_complete.txt
     └── TDHF calculation output
         └── tdhf3d.inp (input configuration)
             └── Initial User Configuration
```

Key Workflow Characteristics:
- Fully automated nuclear structure calculation
- Flexible configuration options
- Comprehensive result tracking
- Comparative data analysis
- Automated reporting