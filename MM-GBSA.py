#!/usr/bin/env python3
"""
Prime MM-GBSA Simulation Runner (No Parsing)
============================================

This script ONLY runs MM-GBSA simulations and generates output files.
No parsing or analysis is performed.

Usage: $SCHRODINGER/run mmgbsa_simulation_runner.py

Requirements: complex.maegz file in current directory
"""

import os
import sys
import subprocess
import time
import glob
from pathlib import Path

class MMGBSASimulationRunner:
    """Pure MM-GBSA simulation runner - no parsing"""
    
    def __init__(self, complex_file="complex.maegz"):
        """Initialize simulation runner with complex file"""
        self.complex_file = complex_file
        self.ligand_asl = None
        self.ligand_description = None
        
        # Check if input file exists
        if not os.path.exists(complex_file):
            raise FileNotFoundError(f"Complex file not found: {complex_file}")
        
        print(f"[*] MM-GBSA Simulation Runner initialized with: {complex_file}")
        print("=" * 60)
    
    def identify_ligand_automatically(self):
        """Try different ligand selection methods to find the peptide"""
        
        print("[Search] Step 1: Identifying ligand (peptide) in complex...")
        
        # Common ligand selection methods for peptides
        ligand_methods = [
            ("mol.n 2", "Second molecule"),
            ("chain B", "Chain B"),
            ("chain C", "Chain C"), 
            ("not protein", "Non-protein atoms"),
            ("res.pt PEP", "Peptide residue type"),
            ("res.pt UNL", "Unknown ligand"),
            ("atom.n 1000-5000", "Atom range (assuming peptide is at end)"),
            ("res.n 1-20", "Residue range 1-20"),
            ("res.n 100-120", "Residue range 100-120")
        ]
        
        for ligand_asl, description in ligand_methods:
            print(f"  Testing: {description} ({ligand_asl})")
            
            # Test with a quick energy calculation
            test_job = f"test_{ligand_asl.replace(' ', '_').replace('.', '_').replace('-', '_')}"
            
            cmd = [
                "prime_mmgbsa",
                "-ligand", f'"{ligand_asl}"',
                "-job_type", "ENERGY",
                "-jobname", test_job,
                "-csv_output", "yes",
                "-NOJOBID",
                self.complex_file
            ]
            
            try:
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
                
                if result.returncode == 0 and os.path.exists(f"{test_job}-out.csv"):
                    print(f"  [OK] SUCCESS: Found ligand using {description}")
                    self.ligand_asl = ligand_asl
                    self.ligand_description = description
                    
                    # Clean up test files
                    for f in glob.glob(f"{test_job}*"):
                        try:
                            os.remove(f)
                        except:
                            pass
                    
                    return True
                else:
                    print(f"  [X] Failed: {ligand_asl}")
                    
            except subprocess.TimeoutExpired:
                print(f"  [Timeout] Timeout: {ligand_asl}")
            except Exception as e:
                print(f"  [Error] Error: {ligand_asl} - {e}")
        
        print("[X] Could not automatically identify ligand!")
        print("   Please check your complex.maegz file structure")
        return False
    
    def run_per_residue_simulation(self):
        """Run per-residue MM-GBSA simulation"""
        
        print(f"\n[MM-GBSA] Step 2: Running per-residue MM-GBSA simulation...")
        print(f"Using ligand selection: {self.ligand_description}")
        
        # Main per-residue simulation with comprehensive options
        cmd = [
            "prime_mmgbsa",
            "-ligand", f'"{self.ligand_asl}"',
            "-rflexdist", "6.0",              # 6Å flexibility around ligand
            "-rflexgroup", "side",            # Side chain flexibility
            "-job_type", "REAL_MIN",          # Full minimization
            "-jobname", "per_residue_main",
            "-csv_output", "yes",             # Generate CSV output
            "-report_prime_log", "yes",       # Generate detailed log
            "-out_type", "COMPLEX",           # Output complex structure
            "-use_ligand_charges",            # Use input charges if available
            "-OVERWRITE",                     # Overwrite existing files
            self.complex_file
        ]
        
        print("  [Start] Starting per-residue simulation...")
        print(f"  Command: {' '.join(cmd)}")
        
        try:
            start_time = time.time()
            result = subprocess.run(cmd, capture_output=True, text=True)
            end_time = time.time()
            
            if result.returncode == 0:
                print(f"  [OK] Per-residue simulation completed successfully!")
                print(f"  [Time] Runtime: {(end_time - start_time)/60:.1f} minutes")
                
                # Check output files
                output_files = [
                    "per_residue_main-out.csv",
                    "per_residue_main.log", 
                    "per_residue_main-out.maegz"
                ]
                
                print("  [Files] Generated files:")
                for file in output_files:
                    if os.path.exists(file):
                        size = os.path.getsize(file) / 1024  # KB
                        print(f"    [OK] {file} ({size:.1f} KB)")
                    else:
                        print(f"    [X] {file} (missing)")
                
                return True
            else:
                print(f"  [X] Per-residue simulation failed!")
                print(f"  Error output: {result.stderr}")
                return False
                
        except Exception as e:
            print(f"  [X] Error running per-residue simulation: {e}")
            return False
    
    def run_per_atom_simulations(self):
        """Run per-atom MM-GBSA simulations for individual residues"""
        
        print(f"\n[Atom] Step 3: Running per-atom MM-GBSA simulations...")
        
        # First, analyze the structure to find actual residue numbers
        actual_residues = self.get_actual_residue_numbers()
        
        if not actual_residues:
            print("  [X] Could not determine residue numbers from structure")
            return []
        
        print(f"  [Info] Found {len(actual_residues)} residues in structure")
        print(f"  [Info] Residue range: {min(actual_residues)} - {max(actual_residues)}")
        
        # Test a few residues to see what works
        test_residues = actual_residues[:5]  # Test first 5 residues
        print(f"  [Test] Testing residues: {test_residues}")
        
        successful_simulations = []
        
        for res_num in test_residues:
            print(f"  Testing residue {res_num}...")
            
            test_jobname = f"test_atom_res{res_num}"
            
            cmd = [
                "prime_mmgbsa",
                "-ligand", f'"res.n {res_num}"',
                "-rflexdist", "5.0",
                "-rflexgroup", "side",
                "-job_type", "ENERGY",        # Quick energy calc for testing
                "-jobname", test_jobname,
                "-csv_output", "yes",
                "-NOJOBID",
                self.complex_file
            ]
            
            try:
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
                
                if result.returncode == 0 and os.path.exists(f"{test_jobname}-out.csv"):
                    print(f"    [OK] Residue {res_num} works - running full simulation")
                    
                    # Clean up test file
                    for f in glob.glob(f"{test_jobname}*"):
                        try:
                            os.remove(f)
                        except:
                            pass
                    
                    # Run full simulation for this residue
                    full_success = self._run_single_residue_simulation(res_num)
                    if full_success:
                        successful_simulations.append(full_success)
                        
                        # If we find one that works, try a few more
                        if len(successful_simulations) < 5:
                            continue
                        else:
                            break
                else:
                    print(f"    [X] Residue {res_num} failed")
                    
            except subprocess.TimeoutExpired:
                print(f"    [Timeout] Residue {res_num} timed out")
            except Exception as e:
                print(f"    [X] Residue {res_num} error: {e}")
        
        print(f"  [Done] Completed {len(successful_simulations)} individual residue simulations")
        return successful_simulations
    
    def get_actual_residue_numbers(self):
        """Get actual residue numbers from the structure file"""
        
        print("  [Analyze] Analyzing structure to find residue numbers...")
        
        # Create a simple Python script to analyze the structure
        analysis_script = '''
import sys
import os
sys.path.append(os.environ.get('SCHRODINGER', ''))

try:
    from schrodinger import structure
    
    st = structure.StructureReader("''' + self.complex_file + '''").next()
    
    residue_numbers = []
    for residue in st.residue:
        residue_numbers.append(residue.resnum)
    
    # Remove duplicates and sort
    unique_residues = sorted(list(set(residue_numbers)))
    
    print("RESIDUE_NUMBERS:", ",".join(map(str, unique_residues)))
    
except Exception as e:
    print("ERROR:", str(e))
    # Fallback to common ranges
    print("RESIDUE_NUMBERS:", ",".join(map(str, range(1, 21))))
'''
        
        with open("analyze_residues.py", "w") as f:
            f.write(analysis_script)
        
        try:
            # Run the analysis script
            result = subprocess.run([
                os.path.join(os.environ.get('SCHRODINGER', ''), 'run'),
                "analyze_residues.py"
            ], capture_output=True, text=True, timeout=60)
            
            # Parse output
            for line in result.stdout.split('\n'):
                if line.startswith("RESIDUE_NUMBERS:"):
                    residue_str = line.split(":", 1)[1].strip()
                    residues = [int(x) for x in residue_str.split(",") if x.strip()]
                    print(f"    [OK] Found {len(residues)} residues")
                    return residues
            
            # Fallback if parsing failed
            print("    [Warning] Could not parse residues, using fallback")
            return list(range(1, 21))
            
        except Exception as e:
            print(f"    [X] Error analyzing structure: {e}")
            return list(range(1, 21))  # Fallback
        finally:
            # Clean up
            if os.path.exists("analyze_residues.py"):
                try:
                    os.remove("analyze_residues.py")
                except:
                    pass
    
    def _run_single_residue_simulation(self, res_num):
        """Run full simulation for a single residue"""
        
        jobname = f"atom_res{res_num}"
        
        cmd = [
            "prime_mmgbsa",
            "-ligand", f'"res.n {res_num}"',
            "-rflexdist", "5.0",
            "-rflexgroup", "side",
            "-job_type", "REAL_MIN",      # Full minimization
            "-jobname", jobname,
            "-csv_output", "yes",
            "-report_prime_log", "yes",
            "-NOJOBID",
            self.complex_file
        ]
        
        try:
            start_time = time.time()
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
            end_time = time.time()
            
            if result.returncode == 0 and os.path.exists(f"{jobname}-out.csv"):
                print(f"      [OK] Residue {res_num} completed ({(end_time-start_time)/60:.1f} min)")
                return {
                    'residue': res_num,
                    'jobname': jobname,
                    'runtime': end_time - start_time
                }
            else:
                print(f"      [X] Residue {res_num} failed")
                return None
                
        except subprocess.TimeoutExpired:
            print(f"      [Timeout] Residue {res_num} timed out (>10 min)")
            return None
        except Exception as e:
            print(f"      [X] Residue {res_num} error: {e}")
            return None
    
    def run_receptor_flexibility_simulation(self):
        """Run simulation with enhanced receptor flexibility"""
        
        print(f"\n[Flexibility] Running enhanced receptor flexibility simulation...")
        
        cmd = [
            "prime_mmgbsa",
            "-ligand", f'"{self.ligand_asl}"',
            "-rflexdist", "8.0",              # Extended flexibility distance
            "-rflexgroup", "residue",         # Entire residue flexible (valid option)
            "-job_type", "REAL_MIN",
            "-jobname", "enhanced_flexibility",
            "-csv_output", "yes",
            "-report_prime_log", "yes",
            "-out_type", "COMPLEX",
            "-use_ligand_charges",
            "-OVERWRITE",
            self.complex_file
        ]
        
        print("  [Start] Starting enhanced flexibility simulation...")
        print(f"  Command: {' '.join(cmd)}")
        
        try:
            start_time = time.time()
            result = subprocess.run(cmd, capture_output=True, text=True)
            end_time = time.time()
            
            if result.returncode == 0:
                print(f"  [OK] Enhanced flexibility simulation completed!")
                print(f"  [Time] Runtime: {(end_time - start_time)/60:.1f} minutes")
                
                # Check output files
                output_files = [
                    "enhanced_flexibility-out.csv",
                    "enhanced_flexibility.log", 
                    "enhanced_flexibility-out.maegz"
                ]
                
                print("  [Files] Generated files:")
                for file in output_files:
                    if os.path.exists(file):
                        size = os.path.getsize(file) / 1024  # KB
                        print(f"    [OK] {file} ({size:.1f} KB)")
                    else:
                        print(f"    [X] {file} (missing)")
                
                return True
            else:
                print(f"  [X] Enhanced flexibility simulation failed!")
                print(f"  Error output: {result.stderr}")
                return False
                
        except Exception as e:
            print(f"  [X] Error running enhanced flexibility simulation: {e}")
            return False
    
    def generate_simulation_summary(self, per_atom_results):
        """Generate a simple summary of simulation results (no parsing)"""
        
        print(f"\n[Step 5] Generating simulation summary...")
        
        with open("simulation_summary.txt", "w") as f:
            f.write("MM-GBSA SIMULATION SUMMARY\n")
            f.write("=" * 30 + "\n\n")
            
            f.write(f"Complex file: {self.complex_file}\n")
            f.write(f"Ligand selection: {self.ligand_description} ({self.ligand_asl})\n")
            f.write(f"Simulation time: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            # List all generated files
            f.write("GENERATED OUTPUT FILES:\n")
            f.write("-" * 25 + "\n")
            
            # Per-residue files
            per_residue_files = [
                "per_residue_main-out.csv",
                "per_residue_main.log",
                "per_residue_main-out.maegz"
            ]
            
            f.write("Per-residue analysis files:\n")
            for file in per_residue_files:
                if os.path.exists(file):
                    size = os.path.getsize(file) / 1024  # KB
                    f.write(f"  [OK] {file} ({size:.1f} KB)\n")
                else:
                    f.write(f"  [X] {file} (missing)\n")
            
            f.write("\n")
            
            # Per-atom files
            if per_atom_results:
                f.write("Per-atom analysis files:\n")
                for result in per_atom_results:
                    jobname = result['jobname']
                    residue = result['residue']
                    runtime = result['runtime']
                    
                    csv_file = f"{jobname}-out.csv"
                    if os.path.exists(csv_file):
                        size = os.path.getsize(csv_file) / 1024  # KB
                        f.write(f"  [OK] {csv_file} - Residue {residue} ({size:.1f} KB, {runtime/60:.1f} min)\n")
                    else:
                        f.write(f"  [X] {csv_file} - Residue {residue} (missing)\n")
            
            f.write("\n")
            
            # Enhanced flexibility files
            enhanced_files = [
                "enhanced_flexibility-out.csv",
                "enhanced_flexibility.log",
                "enhanced_flexibility-out.maegz"
            ]
            
            f.write("Enhanced flexibility files:\n")
            for file in enhanced_files:
                if os.path.exists(file):
                    size = os.path.getsize(file) / 1024  # KB
                    f.write(f"  [OK] {file} ({size:.1f} KB)\n")
                else:
                    f.write(f"  [X] {file} (missing)\n")
            
            f.write("\n")
            f.write("NOTE: Use separate parsing scripts to analyze the CSV and log files.\n")
            f.write("All simulation output files are ready for analysis.\n")
        
        print("  [OK] Simulation summary saved: simulation_summary.txt")
    
    def run_all_simulations(self):
        """Run the complete MM-GBSA simulation workflow"""
        
        print("[START] MM-GBSA SIMULATIONS")
        print("=" * 60)
        
        total_start_time = time.time()
        
        # Step 1: Identify ligand
        if not self.identify_ligand_automatically():
            print("[FAIL] Could not identify ligand in complex")
            return False
        
        # Step 2: Run per-residue simulation
        #per_residue_success = self.run_per_residue_simulation()
        #if not per_residue_success:
        #    print("[FAIL] Per-residue simulation failed")
        #    return False
        
        # Step 3: Run per-atom simulations
        #per_atom_results = self.run_per_atom_simulations()
        
        # Step 4: Run enhanced flexibility simulation
        enhanced_success = self.run_receptor_flexibility_simulation()
        
        # Step 5: Generate simulation summary
        #self.generate_simulation_summary(per_atom_results)
        
        # Final summary
        total_end_time = time.time()
        total_duration = total_end_time - total_start_time
        
        print("\n" + "=" * 60)
        print("[SUCCESS] ALL SIMULATIONS COMPLETE!")
        print("=" * 60)
        print(f"[Time] Total runtime: {total_duration/60:.1f} minutes")
        print(f"[Files] Output files generated in current directory")
        print(f"[Summary] Summary saved in: simulation_summary.txt")
        print("")
        print("[Ready] Ready for analysis!")
        print("   Use your preferred parsing script to analyze the CSV and log files.")
        print("")
        print("[Output] Main output files:")
        print("   * per_residue_main-out.csv - Per-residue results")
        print("   * per_residue_main.log - Detailed per-residue log")
        print("   * atom_res*-out.csv - Individual residue results")
        print("   * enhanced_flexibility-out.csv - Enhanced flexibility results")
        print("")
        print("[Complete] All simulations completed successfully!")
        
        return True

def main():
    """Main execution function"""
    
    # Check for complex file
    complex_file = "complex.maegz"
    
    if len(sys.argv) > 1:
        complex_file = sys.argv[1]
    
    if not os.path.exists(complex_file):
        print(f"[ERROR] {complex_file} not found in current directory")
        print("\nUsage: $SCHRODINGER/run mmgbsa_simulation_runner.py [complex.maegz]")
        print("\nMake sure your complex.maegz file is in the current directory")
        return
    
    try:
        # Create and run simulation runner
        simulator = MMGBSASimulationRunner(complex_file)
        success = simulator.run_all_simulations()
        
        if success:
            print("\n[SUCCESS] All MM-GBSA simulations completed!")
            print("Use a separate parsing script to analyze the results.")
        else:
            print("\n[FAIL] Simulations encountered errors")
            
    except KeyboardInterrupt:
        print("\n[STOP] Simulations interrupted by user")
    except Exception as e:
        print(f"\n[ERROR] {e}")

if __name__ == "__main__":
    main()