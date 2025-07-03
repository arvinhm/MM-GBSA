# MM-GBSA Analysis Pipeline

A comprehensive Python pipeline for running MM-GBSA (Molecular Mechanics Generalized Born Surface Area) simulations and analyzing per-residue energy contributions for protein-peptide complexes.

## 📋 Table of Contents

- [Overview](#overview)
- [Prerequisites](#prerequisites)
- [Installation & Setup](#installation--setup)
- [Quick Start](#quick-start)
- [Detailed Usage](#detailed-usage)
- [Output Files](#output-files)
- [Analysis Results](#analysis-results)
- [Troubleshooting](#troubleshooting)
- [Examples](#examples)

## 🔬 Overview

This pipeline consists of two main components:

1. **`MM-GBSA.py`** - Automated MM-GBSA simulation runner
2. **`parsing.py`** - Log file parser and energy analysis tool

### Key Features

- ✅ **Automatic ligand detection** for protein-peptide complexes
- ✅ **Per-residue energy decomposition** analysis
- ✅ **Automatic peptide vs protein separation** based on residue patterns
- ✅ **Multiple simulation types** (per-residue, per-atom, enhanced flexibility)
- ✅ **Energy minimization** with receptor flexibility
- ✅ **Comprehensive output parsing** with CSV generation
- ✅ **Most negative energy selection** across multiple iterations

## 🛠 Prerequisites

### Software Requirements

1. **Schrödinger Suite** (2023-2 or later)
   - Prime MM-GBSA module
   - Python environment with Schrödinger libraries

2. **Python Packages**
   ```bash
   pip install pandas numpy
   ```

### Input Files

- **`complex.maegz`** - Protein-peptide complex structure in Maestro format

## ⚙️ Installation & Setup

1. **Clone or download the scripts**
   ```bash
   # Place both scripts in your working directory
   MM-GBSA.py
   parsing.py
   ```

2. **Set up Schrödinger environment**
   ```bash
   # Make sure Schrödinger is properly installed and licensed
   echo $SCHRODINGER
   # Should return path like: /opt/schrodinger2023-2
   ```

3. **Prepare your complex structure**
   ```bash
   # Ensure your complex.maegz file is in the working directory
   ls -la complex.maegz
   ```

## 🚀 Quick Start

### Step 1: Run MM-GBSA Simulations

```bash
# Run the complete MM-GBSA simulation pipeline
$SCHRODINGER/run MM-GBSA.py

# Or specify a different complex file
$SCHRODINGER/run MM-GBSA.py my_complex.maegz
```

### Step 2: Parse Results and Generate Analysis

```bash
# Parse the log files and generate CSV analysis
python parsing.py
```

### That's it! 

You'll get two CSV files with the most negative energy residues:
- `peptide_residue_analysis.csv`
- `protein_residue_analysis.csv`

## 📖 Detailed Usage

### MM-GBSA Simulation Runner (`MM-GBSA.py`)

#### What it does:
1. **Automatic Ligand Detection** - Tests multiple ASL expressions to find peptide
2. **Per-Residue Analysis** - Complete energy breakdown for all residues
3. **Per-Atom Analysis** - Individual residue simulations (subset)
4. **Enhanced Flexibility** - Extended receptor flexibility simulation
5. **Summary Generation** - Creates simulation summary report

#### Command Line Options:
```bash
# Basic usage
$SCHRODINGER/run MM-GBSA.py

# Custom complex file
$SCHRODINGER/run MM-GBSA.py path/to/complex.maegz

# The script will automatically:
# - Detect ligand (peptide) in your complex
# - Run multiple simulation types
# - Generate output files
# - Create simulation summary
```

#### Simulation Parameters:
- **Receptor flexibility**: 6-8 Å around ligand
- **Flexibility type**: Side chain or full residue
- **Minimization**: Full energy minimization
- **Solvent model**: Generalized Born (GB)
- **Output format**: CSV + detailed logs

### Energy Analysis Parser (`parsing.py`)

#### What it does:
1. **Parses MM-GBSA logs** - Extracts all energy breakdown tables
2. **Handles multiple iterations** - Finds most negative energy per residue
3. **Separates peptide/protein** - Automatic classification
4. **Generates CSV files** - Clean, analysis-ready data
5. **Provides summary statistics** - Top contributing residues

#### Usage:
```bash
python parsing.py
```

#### Peptide vs Protein Classification:
The script automatically identifies peptide residues based on:
- **Residue types**: ACE (N-terminus), NMA (C-terminus)
- **Residue numbers**: Typically 1-10 for peptides
- **Common peptide residues**: ILE, ASP, CYS, GLN in low-numbered positions

## 📁 Output Files

### MM-GBSA Simulation Outputs

| File | Description | Size |
|------|-------------|------|
| `per_residue_main.log` | Detailed per-residue energy breakdown | ~200 KB |
| `per_residue_main-out.csv` | Per-residue CSV summary | ~3 KB |
| `per_residue_main-out.maegz` | Optimized complex structure | ~600 KB |
| `enhanced_flexibility.log` | Enhanced flexibility analysis | ~60 KB |
| `atom_res*-out.csv` | Individual residue analyses | ~1-2 KB each |
| `simulation_summary.txt` | Complete simulation report | ~2 KB |

### Analysis Outputs

| File | Description | Contents |
|------|-------------|----------|
| `peptide_residue_analysis.csv` | Peptide energy analysis | 7 peptide residues with most negative energies |
| `protein_residue_analysis.csv` | Protein energy analysis | 226 protein residues with most negative energies |

### CSV File Format

Both CSV files contain these columns:

```csv
res_full,res_type,res_num,bond,angle,torsional,el,lj,sgb,nonpolar,lipo,hbond,pipack,selfcont,total
A:GLN_6,GLN,6,0.54,1.81,3.92,-61.28,-3.68,-3.25,0.0,-1.26,-0.59,-0.08,0.0,-63.86
```

**Column Descriptions:**
- `res_full`: Full residue identifier (e.g., A:GLN_6)
- `res_type`: Amino acid type (e.g., GLN)
- `res_num`: Residue number
- `bond/angle/torsional`: Bonded energy terms (kcal/mol)
- `el`: Electrostatic energy (kcal/mol)
- `lj`: Lennard-Jones (van der Waals) energy (kcal/mol)
- `sgb`: Generalized Born solvation energy (kcal/mol)
- `lipo/hbond/pipack`: Additional interaction energies (kcal/mol)
- `total`: Total interaction energy (kcal/mol)

## 📊 Analysis Results

### Example Output Summary

```
Parsing per_residue_main.log...
============================================================
Found data for 233 unique residues
Iterations per residue: 3 - 4
Peptide residues identified: 7
Protein residues identified: 226

Peptide residue energy breakdown:
================================================================================
   A:ACE_1 | ACE |   1 | Total:   -30.43 | El:   -26.86 | LJ:    -1.63 | SGB:    -0.83
   A:GLN_6 | GLN |   6 | Total:   -63.86 | El:   -61.28 | LJ:    -3.68 | SGB:    -3.25

Top 10 most negative protein residues:
================================================================================
   A:GLN_107 | GLN | 107 | Total:   -81.06 | El:   -68.82 | LJ:    -7.52 | SGB:    -3.08
   A:GLN_134 | GLN | 134 | Total:   -77.77 | El:   -71.85 | LJ:    -4.32 | SGB:    -1.44

Summary:
Most negative peptide residue: A:GLN_6 (Total: -63.86)
Most negative protein residue: A:GLN_107 (Total: -81.06)
```

### Interpretation Guide

- **Negative energies** = Favorable interactions
- **Positive energies** = Unfavorable interactions
- **Electrostatic (el)** = Charge-charge interactions
- **Lennard-Jones (lj)** = Van der Waals interactions
- **SGB** = Solvation energy contribution

## 🔧 Troubleshooting

### Common Issues and Solutions

#### Issue: "Complex file not found"
```bash
[ERROR] complex.maegz not found in current directory
```
**Solution**: Ensure your complex file is in the working directory
```bash
ls -la *.maegz
# Should show your complex file
```

#### Issue: "Could not automatically identify ligand"
```bash
[X] Could not automatically identify ligand!
```
**Solutions**:
1. Check that your complex contains both protein and peptide
2. Verify peptide is properly prepared in Maestro
3. Try manual ligand selection by modifying the script

#### Issue: "No residue data found in log file"
```bash
No residue data found in the log file
```
**Solutions**:
1. Check that `per_residue_main.log` exists and is not empty
2. Verify MM-GBSA simulation completed successfully
3. Look for error messages in the log file

#### Issue: Schrödinger environment not found
```bash
bash: prime_mmgbsa: command not found
```
**Solutions**:
1. Source Schrödinger environment:
   ```bash
   source /opt/schrodinger2023-2/schrodinger.bashrc
   ```
2. Check installation:
   ```bash
   echo $SCHRODINGER
   which prime_mmgbsa
   ```

### Performance Tips

1. **Large systems**: Simulations may take 30-60 minutes
2. **Memory usage**: Ensure at least 8 GB RAM available
3. **Disk space**: Reserve 1-2 GB for output files
4. **CPU cores**: Utilizes available cores automatically

## 📝 Examples

### Example 1: Basic Protein-Peptide Analysis

```bash
# 1. Prepare your complex
cp my_protein_peptide_complex.maegz complex.maegz

# 2. Run simulations
$SCHRODINGER/run MM-GBSA.py

# 3. Parse results
python parsing.py

# 4. Analyze CSV files
head peptide_residue_analysis.csv
head protein_residue_analysis.csv
```

### Example 2: Custom Complex File

```bash
# Run with custom file name
$SCHRODINGER/run MM-GBSA.py protein_ligand.maegz

# Parsing will still use per_residue_main.log
python parsing.py
```

### Example 3: Batch Processing

```bash
# Process multiple complexes
for complex in complex_*.maegz; do
    echo "Processing $complex..."
    $SCHRODINGER/run MM-GBSA.py $complex
    python parsing.py
    
    # Rename outputs
    base=$(basename $complex .maegz)
    mv peptide_residue_analysis.csv ${base}_peptide.csv
    mv protein_residue_analysis.csv ${base}_protein.csv
done
```

## 🎯 Expected Results

### Typical Runtime
- **Small systems** (< 100 residues): 15-30 minutes
- **Medium systems** (100-300 residues): 30-60 minutes
- **Large systems** (> 300 residues): 1-2 hours

### File Sizes
- **Log files**: 50-500 KB
- **CSV files**: 1-10 KB
- **Structure files**: 100 KB - 1 MB

### Success Indicators
✅ All simulation steps complete without errors
✅ Log files contain energy breakdown tables
✅ CSV files generated with expected residue counts
✅ Most negative energies identified for key residues

---

## 🆘 Support

If you encounter issues:

1. **Check the simulation summary**: `simulation_summary.txt`
2. **Review log files** for error messages
3. **Verify Schrödinger installation** and licensing
4. **Ensure complex structure quality** (no missing atoms, proper protonation)

For additional help, check the Schrödinger documentation or contact your system administrator.

---

**Happy Analyzing! 🧬💻** 
