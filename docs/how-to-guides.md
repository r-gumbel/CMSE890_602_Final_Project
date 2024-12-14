# How-to Guides

## Configuring the Workflow from the Command Line

### Direct Snakemake Configuration

You can directly configure the workflow using the `--config` flag when running Snakemake:

```bash
# Basic usage
snakemake -j1 --config A=48 Z=20

# Adding more parameters
snakemake -j1 --config A=48 Z=20 email=user@example.com skyrme=SLy4dL
```

#### Configuration Options

| Parameter | Description | Example |
|-----------|-------------|---------|
| `A` | Atomic Mass | 48 |
| `Z` | Atomic Number | 20 |
| `email` | Notification Email | user@example.com |
| `skyrme` | Skyrme Force | SLy4dL |

### Overriding Configuration File

Command-line configurations take precedence over the `config.yaml` file, allowing for quick, flexible parameter changes without modifying the configuration file.

## Changing Nucleus Parameters

To calculate for a different nucleus, use the command-line configuration:

```yaml
# Command line
snakemake -j1 --config A=208 Z=82
```

## Configuring Cluster Resources

Adjust SLURM parameters directly in the configuration:

```bash
snakemake -j1 --config A=40 Z=20 time_limit=48:00:00 mem_per_cpu=8GB
```

## Customizing Module Loads

Add or modify modules using the command line:

```bash
snakemake -j1 --config A=40 Z=20 modules='["intel/2020", "intel-mpi/2020", "python/3.8"]'
```

## How to Send Job Notifications

Specify email directly during workflow execution:

```bash
snakemake -j1 --config A=40 Z=20 email=researcher@university.edu
```

## How to Modify Convergence Parameters

While detailed convergence parameters are typically set in the input file, you can pass basic convergence hints:

```bash
snakemake -j1 --config A=40 Z=20 convergence_tolerance=1e-4
```

## Debugging Common Issues

1. **Job Not Submitting**
   - Verify all required parameters
   - Check module availability
   - Ensure file paths are correct

2. **Long Calculation Times**
   - Reduce grid size
   - Lower convergence precision
   - Increase computational resources

3. **Missing Modules**
   - Contact cluster administrator
   - Use alternative module versions

### Debugging Command-Line Options

```bash
# Dry run to see what would happen
snakemake -j1 --config A=40 Z=20 --dry-run

# Verbose output for detailed information
snakemake -j1 --config A=40 Z=20 -p
```

## Advanced Configuration: Directly Editing tdhf3d.inp

For the most granular control over your nuclear structure calculation, you can directly modify the `tdhf3d.inp` template file.

### Locating the Template

The template is typically found in:
```
VU-TDHF3D.template/run/tdhf3d.inp
```

### Key Configuration Sections

#### Convergence Parameters

```
7500  2500  6.0D-4 7.0D-2   # itrbx,mtrbx,serr,derr
```
- `itrbx`: Maximum iteration count
- `mtrbx`: Maximum iteration box
- `serr`: Convergence error tolerance
- `derr`: Additional error tolerance

#### Nuclei Specification

```
40.0D0 12.0D0   # Mass (fmass)
20.0D0 6.0D0    # Charge (fcharg)
```

#### Computational Domain

```
3.0 3.0   # radinf(1,if) - Radial infinity parameters
3.0 3.0   # radinf(2,if)
3.0 3.0   # radinf(3,if)
```

#### Advanced Numerical Controls

```
0.02   35.0   # x0dmp,e0dmp - Damping parameters
```

### Recommended Workflow

1. Copy the template file:
```bash
cp VU-TDHF3D.template/run/tdhf3d.inp custom_tdhf3d.inp
```

2. Make your modifications

3. Carefully test changes incrementally

### Warning

⚠️ **Caution**: 
- Incorrect modifications can lead to:
  - Calculation failures
  - Unrealistic nuclear models
  - Extremely long computation times

### Best Practices

- Document your changes
- Keep a log of parameter modifications
- Compare results with standard configurations
- Validate against known nuclear structure data

### Example Advanced Configurations

#### Increasing Precision
```
# Increase iteration count and lower error tolerance
10000  3000  1.0D-5 5.0D-3   # More precise calculation
```

#### Exploring Constrained Calculations
```
# Enable constraint parameters
1   20   00  20   # Different print and constraint settings
0 0.000 20  1.90 5e-5   # Constrained calculation parameters
```

### Recommended Reading

- VU-TDHF3D User Manual
- Nuclear structure computational guides
- Peer-reviewed papers using similar methodologies

## Getting Help

- Consult your research advisor
- Reach out to VU-TDHF3D developers
- Compare with established computational nuclear physics workflows