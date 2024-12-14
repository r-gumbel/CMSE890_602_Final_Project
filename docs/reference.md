# Technical Reference

## Configuration Parameters

### Nucleus Configuration

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `A` | Integer | Mass Number | 40 |
| `Z` | Integer | Atomic Number | 20 |

### SLURM Configuration

| Parameter | Type | Description | Example |
|-----------|------|-------------|---------|
| `time_limit` | String | Maximum wall clock time | "24:00:00" |
| `ntasks` | Integer | Number of tasks | 1 |
| `nodes` | Integer | Number of compute nodes | 1 |
| `mem_per_cpu` | String | Memory per CPU | "4GB" |
| `cpus-per-task` | Integer | CPUs per task | 80 |
| `mail_user` | String | Notification email | "user@example.com" |

## Input File Parameters

### TDHF3D Input File

| Parameter | Description | Default/Example |
|-----------|-------------|-----------------|
| `nof` | Number of frames | 1 |
| `imode` | Calculation mode | 0 |
| `serr` | Convergence tolerance | 6.0D-4 |
| `derr` | Error tolerance | 7.0D-2 |
| `itrbx` | Iteration box | 7500 |
| `mtrbx` | Maximum iteration box | 2500 |

## Module Dependencies

### Recommended Modules

- `intel/2020`
- `intel-mpi/2020`
- `python/3.8` (optional)

## File Structure

```
project_root/
│
├── config.yaml
├── snakefiles/
│   ├── Snakefile
│   ├── tdhf.smk
│   └── report.smk
│
└── logs/
    └── A_[A]_Z_[Z]_[Skyrme]/
        ├── report.txt
        ├── tdhf_status.txt
        └── job_complete.txt
```

## Performance Metrics

Tracked in report:
- Total Energy
- Quadrupole Moment (Q20)
- Deformation Parameter (β)
- Relative Errors
- Convergence Iterations

## Error Handling

- Job status tracking
- Email notifications
- Comprehensive logging

## External Dependencies

- Snakemake
- SLURM
- VU-TDHF3D code
- Python 3.x
- NumPy

## Computational Complexity

Scaling factors:
- Nucleus size (A)
- Computational resources
- Convergence tolerance