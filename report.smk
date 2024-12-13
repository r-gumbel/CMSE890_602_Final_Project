# Snakefile
"""
This workflow generates reports from completed VU-TDHF3D 
calculations. It is triggered automatically after the TDHF job
completes.
"""
import subprocess
from pathlib import Path
from datetime import datetime
import numpy as np  # NEW: Required for beta calculation

def beta(A: int | float, Z: int | float, Q: float) -> float:
    """Function to calculate deformation parameter `beta` of a nuclei in its ground state.

    Args:
        A (int or float): atomic mass (i.e. total no. nucleons)
        Z (int or float): atomic number (i.e. total number of protons)
        Q (float): quadrupole moment of nucleus in its ground state.  
    """
    R = 1.2 * A**(1/3)
    return np.sqrt(5 * np.pi) / (3 * Z * R**2) * Q

# [Previous functions remain unchanged: parse_hfb_data, get_energy, get_Q20, parse_output_file]

def calculate_relative_error(tdhf_val, hfb_val):
    """
    Calculate relative error between TDHF and HFB values
    relative_error = |experimental - theoretical|/|theoretical| * 100%
    """
    try:
        # Convert strings to floats if needed
        tdhf_num = float(tdhf_val) if isinstance(tdhf_val, str) else tdhf_val
        hfb_num = float(hfb_val) if isinstance(hfb_val, str) else hfb_val
        
        if hfb_num == 0:  # Avoid division by zero
            return None
            
        rel_error = abs(tdhf_num - hfb_num) / abs(hfb_num) * 100.0
        return rel_error
        
    except (ValueError, TypeError) as e:
        print(f"Error calculating relative error: {e}")
        return None

rule generate_report:
    input:
        job_complete = f"logs/{run_id}/job_complete.txt"
    output:
        report = f"logs/{run_id}/report.txt"
    run:
        print("Generating report...")
        
        # Get the numerical values
        tdhf_energy, tdhf_q20 = parse_output_file(working_dir, run_dir, config)
        hfb_data = parse_hfb_data(config['nucleus']['A'], config['nucleus']['Z'])
        
        # Calculate beta from Q20
        tdhf_beta = None
        if tdhf_q20 is not None:
            tdhf_beta = beta(config['nucleus']['A'], config['nucleus']['Z'], tdhf_q20)
        
        # Calculate relative errors
        energy_error = None
        deformation_error = None
        
        if hfb_data and tdhf_energy and tdhf_beta:
            energy_error = calculate_relative_error(tdhf_energy, hfb_data['energy'])
            deformation_error = calculate_relative_error(tdhf_beta, hfb_data['beta'])
        
        # Write report
        with open(output.report, 'w') as f:
            f.write(f"""TDHF Analysis Report
==============================================
Generated: {datetime.now()}
Run ID: {run_id}
Nucleus: A={config['nucleus']['A']}, Z={config['nucleus']['Z']}
Skyrme: {config['skyrme']}

Analysis Results:
-----------------
TDHF Results:
Energy: {tdhf_energy} MeV
Q20: {tdhf_q20}
Beta (converted): {tdhf_beta:.6f} if tdhf_beta else "Not calculated"}

HFB Comparison:
""")
            if hfb_data:
                f.write(f"""Energy: {hfb_data['energy']} MeV
Beta: {hfb_data['beta']}

Benchmarking Results:
--------------------
Energy Relative Error: {energy_error:.2f}% if energy_error else "Not calculated"}
Beta Relative Error: {deformation_error:.2f}% if deformation_error else "Not calculated"}
""")
            else:
                f.write("No matching HFB data found for this nucleus\n")
            
        print("Report generation complete!")