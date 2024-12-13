# Snakefile
"""
This workflow generates reports from completed VU-TDHF3D 
calculations. It is triggered automatically after the TDHF job
completes.
"""
import subprocess
from pathlib import Path
from datetime import datetime
import numpy as np

def beta(A: int | float, Z: int | float, Q: float) -> float:
    """Function to calculate deformation parameter `beta` of a nuclei in its ground state.

    Args:
        A (int or float): atomic mass (i.e. total no. nucleons)
        Z (int or float): atomic number (i.e. total number of protons)
        Q (float): quadrupole moment of nucleus in its ground state.  
    """
    R = 1.2 * A**(1/3)
    return np.sqrt(5 * np.pi) / (3 * Z * R**2) * Q

def parse_hfb_data(A, Z):
    """Parse HFB data file and return matching nucleus data
    
    Args:
        A (int): Mass number
        Z (int): Atomic number
        
    Returns:
        dict: Nucleus data if found, None otherwise
    """
    try:
        with open("HFB.data", 'r') as f:
            # Skip any header rows by checking if the first field is numeric
            for line in f:
                data = line.strip().split()
                if len(data) < 4:  # Skip lines that don't have enough fields
                    continue
                    
                try:
                    current_A = int(data[0])
                    current_Z = int(data[1])
                    
                    if current_A == A and current_Z == Z:
                        return {
                            'A': current_A,
                            'Z': current_Z,
                            'energy': float(data[2]),
                            'beta': float(data[3])
                        }
                except ValueError:
                    # Skip lines that don't have proper numeric values
                    continue
                    
        print(f"No matching data found for A={A}, Z={Z}")
        return None
        
    except FileNotFoundError:
        print("Warning: HFB.data file not found")
        return None
    except Exception as e:
        print(f"Error parsing HFB data: {e}")
        return None

def get_energy(working_dir, run_dir, config):
    """Extract energy value from output file
    
    Args:
        working_dir (str): Working directory path
        run_dir (str): Run directory name
        config (dict): Configuration dictionary
        
    Returns:
        float: Energy value if found, None otherwise
    """
    search_dir = Path(working_dir) / run_dir / "run"
    matching_files = list(search_dir.glob(f"AL_{config['nucleus']['A']}*/tdhf3d.out"))
    
    if not matching_files:
        raise FileNotFoundError(f"No output files found in {search_dir}")
    
    latest_file = max(matching_files, key=lambda p: p.stat().st_mtime)
    
    # Updated grep command to match exact format
    cmd = f"grep 'energy      (mev) =    ' {latest_file}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.stdout:
        # Extract everything after the equals sign and convert to float
        energy_str = result.stdout.strip().split('=')[-1]
        try:
            energy_value = float(energy_str)
            return energy_value
        except (ValueError, IndexError) as e:
            print(f"Warning: Could not parse energy value from: {energy_str}")
            print(f"Error details: {e}")
            # Optionally print the full line for debugging
            print(f"Full matched line: {result.stdout.strip()}")
            return None
    else:
        print(f"Warning: No energy value found in {latest_file}")
        return None

def parse_hfb_data(A, Z):
    """Parse HFB data file and return matching nucleus data
    
    Args:
        A (int): Mass number
        Z (int): Atomic number
        
    Returns:
        dict: Nucleus data if found, None otherwise
    """
    try:
        with open("HFB.data", 'r') as f:
            # Skip any header rows by checking if the first field is numeric
            for line in f:
                data = line.strip().split()
                if len(data) < 4:  # Skip lines that don't have enough fields
                    continue
                    
                try:
                    current_A = int(data[0])
                    current_Z = int(data[1])
                    
                    if current_A == A and current_Z == Z:
                        return {
                            'A': current_A,
                            'Z': current_Z,
                            'energy': float(data[2]),
                            'beta': float(data[3])
                        }
                except ValueError:
                    # Skip lines that don't have proper numeric values
                    continue
                    
        print(f"No matching data found for A={A}, Z={Z}")
        return None
        
    except FileNotFoundError:
        print("Warning: HFB.data file not found")
        return None
    except Exception as e:
        print(f"Error parsing HFB data: {e}")
        return None

def get_Q20(working_dir, run_dir, config):
    """Extract Q20 value from output file"""
    search_dir = Path(working_dir) / run_dir / "run"
    matching_files = list(search_dir.glob(f"AL_{config['nucleus']['A']}*/tdhf3d.out"))
    
    if not matching_files:
        raise FileNotFoundError(f"No output files found in {search_dir}")
    
    latest_file = max(matching_files, key=lambda p: p.stat().st_mtime)
    cmd = f"grep 'total:' {latest_file} | tail -1 | awk '{{print $3}}'"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.stdout:
        try:
            return float(result.stdout.strip())
        except ValueError:
            print(f"Warning: Could not convert Q20 value to float: {result.stdout.strip()}")
            return None
    return None

def parse_output_file(working_dir, run_dir, config):
    """Parse the output file and return energy and Q20 values"""
    try:
        energy = get_energy(working_dir, run_dir, config)
        Q20 = get_Q20(working_dir, run_dir, config)
        return energy, Q20
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return None, None

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
        
        if hfb_data and tdhf_energy and tdhf_beta is not None:
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
Beta (converted): {f'{tdhf_beta:.6f}' if tdhf_beta is not None else 'Not calculated'}

HFB Comparison:
""")
            if hfb_data:
                f.write(f"""Energy: {hfb_data['energy']} MeV
Beta: {hfb_data['beta']}

Benchmarking Results:
--------------------
Energy Relative Error: {f'{energy_error:.2f}%' if energy_error is not None else 'Not calculated'}
Beta Relative Error: {f'{deformation_error:.2f}%' if deformation_error is not None else 'Not calculated'}
""")
            else:
                f.write("No matching HFB data found for this nucleus\n")
            
        print("Report generation complete!")