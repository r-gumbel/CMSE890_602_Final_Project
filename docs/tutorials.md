# Getting Started with Nuclear Structure Calculations

## Introduction

This tutorial will guide you through running your first Time-Dependent Hartree-Fock (TDHF) nuclear structure calculation using our automated workflow script `snake.exe`.

## Prerequisites

Before you begin, ensure you have:
- Snakemake installed
- Access to a SLURM cluster
- The VU-TDHF3D code installed
- The `snake.exe` script in your project directory

## Step 1: Prepare Your Environment

Ensure you are in the project root directory that contains:
- The `Snakefiles` directory
- The `snake.exe` script
- A configuration file (`config.yaml`)

## Step 2: Make the Script Executable

Before first use, make the script executable:

```bash
chmod +x snake.exe
```

## Step 3: Run the Workflow

Execute the script and follow the interactive prompts:

```bash
./snake.exe
```

You will be asked to provide:
1. **Atomic Mass (A)**: Total number of nucleons 
   - Example: 40 for Calcium-40
2. **Atomic Number (Z)**: Number of protons
   - Example: 20 for Calcium
3. **Email Address**: For job completion notifications

### Example Interaction

```
Enter atomic mass A: 40
Enter atomic number Z: 20
Enter email to be notified of job completion: your.email@example.com
Verify email address: your.email@example.com
```

## Step 4: Validating the Output

### Expected Output Directories and Files

#### Simulation Results Directory
The workflow creates a new directory in the parent folder, labeled with:
- Atomic Mass (A)
- Atomic Number (Z)
- Skyrme Parametrization

**Example**: For Calcium-48, you'll see:
`~/A_48_Z_20_SLy4dL`

#### Calculation Output Files
The key output file is `tdhf3d.out`, located in a nested directory:
`A_48_Z_20_SLy4dL/run/AL_48_ZL_20_AR_0_ZR_0_SLy4dL_Ecm0.000_beta0_L0`

#### Report File
A comprehensive summary report is generated in:
`CMSE890_602_Final_Project/logs/A_48_Z_20_SLy4dL/report.txt`

### Understanding the Report

The `report.txt` provides a detailed analysis of the nuclear structure calculation:

```
TDHF Analysis Report
==============================================
Generated: 2024-12-13 15:38:45.412186
Run ID: A_40_Z_20_SLy4dL

Initial Configuration:
---------------------
Nucleus:
- Atomic Mass (A): 40
- Atomic Number (Z): 20
Skyrme Force: SLy4dL
Convergence Tolerance: 0.0006

Analysis Results:
-----------------
TDHF Results:
Energy: -339.6245 MeV
Q20: 1.3667e-05
Beta (converted): 0.000000
Iterations to Convergence: 7500

HFB Comparison:
Energy: -342.689 MeV
Beta: 0.0

Benchmarking Results:
--------------------
Energy Relative Error: 0.89%
Beta Relative Error: Not calculated
```

#### Key Metrics Explained

1. **Energy**: Total energy of the nuclear system in MeV
2. **Q20**: Quadrupole moment
3. **Beta**: Nuclear deformation parameter
4. **Iterations**: Number of computational iterations to converge
5. **Relative Error**: Comparison with Hartree-Fock-Bogoliubov (HFB) model

## Troubleshooting

- Ensure all paths are correct
- Verify email address format
- Check that required modules are available on your cluster
- Confirm the VU-TDHF3D code is properly installed

## Next Steps

- Experiment with different nuclei by rerunning the script
- Try calculations for various atomic masses and numbers
- Compare results with experimental data
- Analyze the convergence and error metrics