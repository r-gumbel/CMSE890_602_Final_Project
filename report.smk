# Snakefile
"""
This workflow generates reports from completed VU-TDHF3D 
calculations. It is triggered automatically after the TDHF job
completes.
"""
import subprocess
from pathlib import Path
from datetime import datetime

def get_energy(working_dir, run_dir, config):
    """Extract energy value from output file"""
    # Construct the path pattern more safely
    search_dir = Path(working_dir) / run_dir / "run"
    # Use Path.glob to find all matching files
    matching_files = list(search_dir.glob(f"AL_{config['nucleus']['A']}*/tdhf3d.out"))
    
    if not matching_files:
        raise FileNotFoundError(f"No output files found in {search_dir}")
    
    # Get the latest file
    latest_file = max(matching_files, key=lambda p: p.stat().st_mtime)
    
    # Use grep through subprocess with the actual file path
    cmd = f"grep 'energy (mev)' {latest_file}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    # Get the last line of the output
    if result.stdout:
        energy = result.stdout.strip().split('\n')[-1]
        return energy
    return None

def get_Q20(working_dir, run_dir, config):
    """Extract Q20 value from output file"""
    # Construct the path pattern more safely
    search_dir = Path(working_dir) / run_dir / "run"
    # Use Path.glob to find all matching files
    matching_files = list(search_dir.glob(f"AL_{config['nucleus']['A']}*/tdhf3d.out"))
    
    if not matching_files:
        raise FileNotFoundError(f"No output files found in {search_dir}")
    
    # Get the latest file
    latest_file = max(matching_files, key=lambda p: p.stat().st_mtime)
    
    # Use grep and awk through subprocess with the actual file path
    cmd = f"grep 'total:' {latest_file} | tail -1 | awk '{{print $3}}'"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    # Return the stripped output
    if result.stdout:
        return result.stdout.strip()
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

def debug_file_search(working_dir, run_dir, config):
    """Debug function to print all found files"""
    search_dir = Path(working_dir) / run_dir / "run"
    print(f"Searching in: {search_dir}")
    print(f"Looking for pattern: AL_{config['nucleus']['A']}*/tdhf3d.out")
    
    all_files = list(search_dir.glob(f"AL_{config['nucleus']['A']}*/tdhf3d.out"))
    print(f"Found {len(all_files)} matching files:")
    for f in all_files:
        print(f"  {f}")

rule generate_report:
    input:
        job_complete = f"logs/{run_id}/job_complete.txt"
    output:
        report = f"logs/{run_id}/report.txt"
    run:
        print("Generating report...")
        
        # Debug file searching (can be removed once everything is working)
        debug_file_search(working_dir, run_dir, config)
        
        # Get the analysis results
        energy, Q20 = parse_output_file(working_dir, run_dir, config)
        
        with open(output.report, 'w') as f:
            f.write(f"""TDHF Analysis Report
==============================================
Generated: {datetime.now()}
Run ID: {run_id}
Nucleus: A={config['nucleus']['A']}, Z={config['nucleus']['Z']}
Skyrme: {config['skyrme']}

Analysis Results:
-----------------
Energy: {energy}
Q20: {Q20}
            """)
        print("Report generation complete")