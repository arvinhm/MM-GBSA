#!/usr/bin/env python3
import re
import pandas as pd
import os
from collections import defaultdict

LOGFILE = "enhanced_flexibility.log"
PEPTIDE_CSV = "peptide_per_residue.csv"
PROTEIN_CSV = "protein_per_residue.csv"

def parse_mmgbsa_log(log_file_path):
    """
    Parse MM-GBSA log file and extract energy breakdown data for residues
    """
    
    # Dictionary to store all energy data for each residue
    # Structure: {residue_id: [list of energy dictionaries]}
    residue_data = defaultdict(list)
    
    # Read the log file
    with open(log_file_path, 'r') as file:
        content = file.read()
    
    # Pattern to match energy breakdown tables
    # Look for the header line followed by residue data
    table_pattern = r'residue\s+bond\s+angle\s+torsional\s+el\s+lj\s+sgb\s+nonpolar\s+lipo\s+hbond\s+pipack\s+selfcont\s+tot\n(.*?)(?=\n\s*\n|\n\s*NOW REPORTING|$)'
    
    tables = re.findall(table_pattern, content, re.DOTALL)
    
    for table in tables:
        lines = table.strip().split('\n')
        
        for line in lines:
            line = line.strip()
            if not line or line.startswith('A:') == False:
                continue
                
            # Parse residue line
            # Format: A:RES_NUM    bond  angle  torsional  el  lj  sgb  nonpolar  lipo  hbond  pipack  selfcont  total
            parts = line.split()
            
            if len(parts) >= 13:  # Ensure we have all required columns
                try:
                    residue_full = parts[0]  # e.g., "A:ILE_16"
                    
                    # Extract residue info
                    residue_match = re.match(r'A:([A-Z]+)_?(\d+)', residue_full)
                    if residue_match:
                        res_type = residue_match.group(1)
                        res_num = residue_match.group(2)
                    else:
                        # Handle cases like A:ACE_1, A:NMA_7
                        res_parts = residue_full.replace('A:', '').split('_')
                        if len(res_parts) == 2:
                            res_type = res_parts[0]
                            res_num = res_parts[1]
                        else:
                            continue
                    
                    # Parse energy values
                    bond = float(parts[1])
                    angle = float(parts[2])
                    torsional = float(parts[3])
                    el = float(parts[4])
                    lj = float(parts[5])
                    sgb = float(parts[6])
                    nonpolar = float(parts[7])
                    lipo = float(parts[8])
                    hbond = float(parts[9])
                    pipack = float(parts[10])
                    selfcont = float(parts[11])
                    total = float(parts[12])
                    
                    # Store the data
                    energy_dict = {
                        'res_full': residue_full,
                        'res_type': res_type,
                        'res_num': int(res_num),
                        'bond': bond,
                        'angle': angle,
                        'torsional': torsional,
                        'el': el,
                        'lj': lj,
                        'sgb': sgb,
                        'nonpolar': nonpolar,
                        'lipo': lipo,
                        'hbond': hbond,
                        'pipack': pipack,
                        'selfcont': selfcont,
                        'total': total
                    }
                    
                    residue_data[residue_full].append(energy_dict)
                    
                except (ValueError, IndexError) as e:
                    print(f"Error parsing line: {line}")
                    print(f"Error: {e}")
                    continue
    
    return residue_data

def get_most_negative_energies(residue_data):
    """
    For each residue, get the entry with the most negative total energy
    """
    most_negative = {}
    
    for residue_id, energy_list in residue_data.items():
        if energy_list:
            # Find the entry with the most negative (minimum) total energy
            min_entry = min(energy_list, key=lambda x: x['total'])
            most_negative[residue_id] = min_entry
    
    return most_negative

def separate_peptide_protein(most_negative_data):
    """
    Separate peptide and protein residues based on residue numbers and types
    """
    peptide_residues = {}
    protein_residues = {}
    
    # Peptide residues are typically: ACE, short peptide sequence, NMA
    # Usually have lower residue numbers (1-10) and include ACE, NMA
    # Protein residues have higher numbers (16+)
    
    for residue_id, data in most_negative_data.items():
        res_type = data['res_type']
        res_num = data['res_num']
        
        # Identify peptide vs protein based on residue type and number
        if (res_type in ['ACE', 'NMA'] or 
            (res_num <= 10 and res_type in ['ILE', 'ASP', 'CYS', 'GLN'])):
            peptide_residues[residue_id] = data
        else:
            protein_residues[residue_id] = data
    
    return peptide_residues, protein_residues

def create_csv_files(peptide_data, protein_data):
    """
    Create CSV files for peptide and protein residues
    """
    
    # Define column order
    columns = ['res_full', 'res_type', 'res_num', 'bond', 'angle', 'torsional', 
               'el', 'lj', 'sgb', 'nonpolar', 'lipo', 'hbond', 'pipack', 'selfcont', 'total']
    
    # Create DataFrames
    if peptide_data:
        peptide_df = pd.DataFrame.from_dict(peptide_data, orient='index')
        peptide_df = peptide_df[columns]
        peptide_df = peptide_df.sort_values('res_num')
        peptide_df.to_csv('peptide_residue_analysis.csv', index=False)
        print(f"Created peptide_residue_analysis.csv with {len(peptide_df)} residues")
        print(f"Peptide residues: {sorted(peptide_df['res_full'].tolist())}")
        
        # Display peptide data
        print("\nPeptide residue energy breakdown:")
        print("="*80)
        for _, row in peptide_df.iterrows():
            print(f"{row['res_full']:>10} | {row['res_type']:>3} | {row['res_num']:>3} | "
                  f"Total: {row['total']:>8.2f} | El: {row['el']:>8.2f} | "
                  f"LJ: {row['lj']:>8.2f} | SGB: {row['sgb']:>8.2f}")
    
    if protein_data:
        protein_df = pd.DataFrame.from_dict(protein_data, orient='index')
        protein_df = protein_df[columns]
        protein_df = protein_df.sort_values('res_num')
        protein_df.to_csv('protein_residue_analysis.csv', index=False)
        print(f"\nCreated protein_residue_analysis.csv with {len(protein_df)} residues")
        print(f"Protein residue range: {min(protein_df['res_num'])} - {max(protein_df['res_num'])}")
        
        # Show top 10 most negative protein residues
        print("\nTop 10 most negative protein residues:")
        print("="*80)
        top_negative = protein_df.nsmallest(10, 'total')
        for _, row in top_negative.iterrows():
            print(f"{row['res_full']:>12} | {row['res_type']:>3} | {row['res_num']:>3} | "
                  f"Total: {row['total']:>8.2f} | El: {row['el']:>8.2f} | "
                  f"LJ: {row['lj']:>8.2f} | SGB: {row['sgb']:>8.2f}")

def main():
    """
    Main function to parse MM-GBSA log and create CSV files
    """
    log_file = "per_residue_main.log"
    
    # Check if log file exists
    if not os.path.exists(log_file):
        print(f"Error: {log_file} not found in current directory")
        return
    
    print(f"Parsing {log_file}...")
    print("="*60)
    
    # Parse the log file
    residue_data = parse_mmgbsa_log(log_file)
    
    if not residue_data:
        print("No residue data found in the log file")
        return
    
    print(f"Found data for {len(residue_data)} unique residues")
    
    # Show how many iterations per residue
    iterations_info = {}
    for residue_id, energy_list in residue_data.items():
        iterations_info[residue_id] = len(energy_list)
    
    max_iterations = max(iterations_info.values())
    min_iterations = min(iterations_info.values())
    print(f"Iterations per residue: {min_iterations} - {max_iterations}")
    
    # Get most negative energies for each residue
    most_negative = get_most_negative_energies(residue_data)
    
    # Separate peptide and protein residues
    peptide_residues, protein_residues = separate_peptide_protein(most_negative)
    
    print(f"Peptide residues identified: {len(peptide_residues)}")
    print(f"Protein residues identified: {len(protein_residues)}")
    
    # Create CSV files
    create_csv_files(peptide_residues, protein_residues)
    
    print("\n" + "="*60)
    print("Analysis complete! Generated files:")
    print("- peptide_residue_analysis.csv")
    print("- protein_residue_analysis.csv")
    
    if peptide_residues and protein_residues:
        # Summary statistics
        peptide_df = pd.DataFrame.from_dict(peptide_residues, orient='index')
        protein_df = pd.DataFrame.from_dict(protein_residues, orient='index')
        
        print(f"\nSummary:")
        print(f"Most negative peptide residue: {peptide_df.loc[peptide_df['total'].idxmin(), 'res_full']} "
              f"(Total: {peptide_df['total'].min():.2f})")
        print(f"Most negative protein residue: {protein_df.loc[protein_df['total'].idxmin(), 'res_full']} "
              f"(Total: {protein_df['total'].min():.2f})")

if __name__ == "__main__":
    main()
