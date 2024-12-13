# Nuclear Structure Calculation Workflow Pseudocode

## Main Workflow (Snakefile)
```plaintext
1. LOAD configuration from config.yaml
2. IF command line arguments provided (A, Z):
   - UPDATE configuration with command line values
3. SET run_id based on nucleus parameters and Skyrme force
4. CREATE log directory structure
5. INCLUDE tdhf.smk and report.smk workflows
6. DEFINE job monitoring functions:
   a. check_job_status:
      - CHECK if job exists in queue
      - RETURN true if job complete
   b. wait_for_job:
      - MONITOR job status every minute
      - RETURN when job completes
7. SET target rule to require final report
8. RUN workflow
```

## TDHF Workflow (tdhf.smk)
```plaintext
1. SET working directory and run directory paths
2. CREATE job name from nucleus parameters
3. DEFINE slurm script generation rule:
   a. GENERATE module load commands
   b. WRITE SLURM script header:
      - Set time limit
      - Set resource requirements
      - Set job name
      - Set email notifications
   c. WRITE environment setup commands:
      - Set OpenMP parameters
      - Set system limits
   d. COPY template directory to run directory
   e. WRITE input file with:
      - Nucleus parameters
      - Calculation settings
      - System configuration
   f. ADD run commands:
      - Execute calculation
      - Record job information
```

## Report Generation Workflow (report.smk)
```plaintext
1. DEFINE output file parsing functions:
   a. get_energy:
      - SEARCH for output file in run directory
      - FIND latest matching file
      - EXTRACT energy value using grep
      - RETURN energy value or None
   
   b. get_Q20:
      - SEARCH for output file in run directory
      - FIND latest matching file
      - EXTRACT Q20 value using grep and awk
      - RETURN Q20 value or None
   
   c. parse_output_file:
      - CALL get_energy
      - CALL get_Q20
      - RETURN both values as tuple
      - HANDLE any file not found errors

2. DEFINE report generation rule:
   a. REQUIRE job completion marker
   b. PARSE output files for results
   c. GENERATE report containing:
      - Timestamp
      - Run information
      - Nucleus parameters
      - Energy results
      - Q20 results
```

## Complete Workflow Execution Order
```plaintext
1. User Input Phase:
   a. READ config.yaml
   b. PROCESS command line arguments
   c. SET final configuration

2. Preparation Phase:
   a. CREATE directory structure
   b. GENERATE run identifiers
   c. PREPARE log directories

3. TDHF Calculation Phase:
   a. GENERATE SLURM script
   b. SUBMIT job to queue
   c. MONITOR job status
   d. WAIT for completion
   e. CREATE completion marker

4. Analysis Phase:
   a. VERIFY job completion
   b. LOCATE output files
   c. EXTRACT relevant data
   d. GENERATE analysis report

5. Completion Phase:
   a. ENSURE all files generated
   b. CLEAN UP temporary files
   c. REPORT final status
```

## Dependency Chain
```plaintext
report.txt
  └── job_complete.txt
       └── TDHF job completion
            └── test_tdhf_{run_id}.slurm
                 └── Configuration
```

This workflow provides:
- Automated job submission and monitoring
- Structured input file generation
- Systematic output analysis
- Comprehensive result reporting

The process runs automatically from submission to analysis, requiring only initial configuration input from the user.