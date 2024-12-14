# Understanding the Nuclear Structure Calculation Workflow

## Workflow Architecture

The Snakemake workflow is designed to automate Time-Dependent Hartree-Fock (TDHF) nuclear structure calculations across multiple computational stages.

### Key Components

1. **Main Snakefile**
   - Coordinates overall workflow
   - Manages configuration
   - Handles job submission and tracking

2. **TDHF Script Generation (tdhf.smk)**
   - Creates SLURM job scripts
   - Configures computational parameters
   - Sets up environment for VU-TDHF3D code

3. **Reporting System (report.smk)**
   - Parses calculation outputs
   - Generates comprehensive analysis reports
   - Compares theoretical results with reference data

## Computational Workflow

```
Config Input → SLURM Script Generation → Job Submission → Calculation → Result Parsing → Report Generation
```

### Computational Parameters

- **Nucleus Specification**
  - Mass Number (A)
  - Atomic Number (Z)

- **Skyrme Force**
  - Defines nuclear interaction model
  - Default: SLy4dL

- **Computational Resources**
  - SLURM cluster configuration
  - Module dependencies
  - Memory and time allocations

## Mathematical Background

### Nuclear Deformation Calculation

The workflow calculates nuclear deformation using:
- Quadrupole moment (Q20)
- Deformation parameter (β)

#### Deformation Parameter Calculation

β = √(5π) / (3Z * R²) * Q20

Where:
- R = 1.2 * A^(1/3)
- Q20 is the quadrupole moment
- Z is the atomic number

## Error Analysis

The report includes:
- Relative errors between theoretical and reference values
- Convergence metrics
- Iteration counts

## Computational Techniques

- Hartree-Fock method
- Time-dependent nuclear structure modeling
- Iterative convergence algorithms