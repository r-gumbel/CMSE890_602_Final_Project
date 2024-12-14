# Nuclear Structure Calculation Workflow

## Overview
A comprehensive workflow for Time-Dependent Hartree-Fock (TDHF) nuclear structure calculations, designed to streamline complex computational physics research.

*"Was Something Left Undone?"*

## Project Structure
```
.
├── Snakefiles/           # Workflow management scripts
├── VU-TDHF3D.template/   # Calculation template
├── snake.exe             # Workflow launcher
├── config.yaml           # Configuration settings
├── HFB.data              # Reference nuclear data
├── logs/                 # Calculation logs
├── outputs/              # Calculation results
└── slurm_scripts/        # Generated SLURM scripts
```

## Motivation
Nuclear interaction modeling presents significant challenges:
- Complex interplay of fundamental forces
- Computationally intensive
- Difficult to calibrate experimentally

## Computational Approach
- **Method**: Time-Dependent Hartree-Fock (TDHF)
- **Tool**: VU-TDHF3D
- **Focus**: Ground state nuclear structure calculations

## Quick Start

1. Configure `config.yaml`
2. Make launcher executable
   ```bash
   chmod +x snake.exe
   ```

3. Run workflow
   ```bash
   ./snake.exe
   ```

## Key Features
- Automated SLURM job generation
- Comprehensive result reporting
- Benchmarking against reference data
- Flexible configuration options

## Usage
Specify nucleus parameters interactively or via command line:
```bash
snakemake -j1 --config A=48 Z=20
```

## Documentation
Detailed guides available in `docs/` directory using MkDocs.