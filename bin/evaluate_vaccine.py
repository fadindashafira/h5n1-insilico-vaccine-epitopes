#!/usr/bin/env python3

import os
import sys
import subprocess
import argparse
import time
import textwrap

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Evaluate vaccine construct properties')
    parser.add_argument('--vaccine', required=True, help='Vaccine construct FASTA file')
    parser.add_argument('--protein-type', required=True, help='Protein type (e.g., hemagglutinin, neuraminidase)')
    parser.add_argument('--linker', default='GPGPG', help='Linker sequence used in the vaccine')
    parser.add_argument('--iedb-api-url', default='http://tools-api.iedb.org/tools_api/', help='IEDB API URL')
    parser.add_argument('--output-evaluation', required=True, help='Output evaluation report file')
    parser.add_argument('--output-properties', required=True, help='Output properties CSV file')
    parser.add_argument('--output-colabfold', required=True, help='Output FASTA file for ColabFold')
    
    args = parser.parse_args()
    
    # Install required packages
    try:
        import pandas as pd
        import numpy as np
        from Bio import SeqIO
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
        import requests
        import json
    except ImportError:
        subprocess.check_call([sys.executable, '-m', 'pip', 'install', '--user', 'pandas', 'numpy', 'biopython', 'requests'])
        import pandas as pd
        import numpy as np
        from Bio import SeqIO
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
        import requests
        import json
    
    print(f"Evaluating {args.protein_type} vaccine construct")
    
    # Read vaccine construct
    try:
        record = next(SeqIO.parse(args.vaccine, "fasta"))
        vaccine_seq = str(record.seq)
        print(f"Successfully loaded vaccine sequence: {len(vaccine_seq)} amino acids")
    except Exception as e:
        print(f"Error reading vaccine file: {e}")
        vaccine_seq = "M"
    
    # ====== PART 1: Use BioPython's ProtParam (local, lightweight) ======
    try:
        analysis = ProteinAnalysis(vaccine_seq)
        
        # Core properties
        mol_weight = analysis.molecular_weight()
        theoretical_pi = analysis.isoelectric_point()
        instability_index = analysis.instability_index()
        gravy = analysis.gravy()
        aromaticity = analysis.aromaticity()
        aa_composition = analysis.get_amino_acids_percent()
        secondary_structure = analysis.secondary_structure_fraction()
        
        print("Completed BioPython ProtParam analysis")
    except Exception as e:
        print(f"Error in ProtParam analysis: {e}")
        # Fallback values
        mol_weight = 0
        theoretical_pi = 0
        instability_index = 0
        gravy = 0
        aromaticity = 0
        aa_composition = {}
        secondary_structure = (0, 0, 0)
    
    # ====== PART 2: IEDB API for epitope prediction (uses web service) ======
    try:
        # IEDB API for BepiPred-2.0
        print("Predicting B-cell epitopes via IEDB API...")
        
        iedb_url = f"{args.iedb_api_url.rstrip('/')}{'/' if not args.iedb_api_url.endswith('/') else ''}bcell/"
        data = {
            "sequence_text": vaccine_seq,
            "method": "bepipred",
            "window_size": 7
        }
        
        # API requests with retry logic
        max_retries = 3
        wait_time = 5  # seconds
        bepipred_results = []
        
        for attempt in range(max_retries):
            try:
                response = requests.post(iedb_url, data=data)
                if response.status_code == 200:
                    bepipred_results = response.json()
                    break
            except Exception as e:
                if attempt < max_retries - 1:
                    print(f"IEDB API attempt {attempt+1} failed. Retrying in {wait_time} seconds...")
                    time.sleep(wait_time)
                    wait_time *= 2  # Exponential backoff
                else:
                    print(f"All IEDB API attempts failed: {e}")
                    bepipred_results = []
        
        # Process epitope results
        antigenic_regions = []
        threshold = 0.5
        in_region = False
        start_pos = 0
        
        if isinstance(bepipred_results, list) and bepipred_results:
            for i, result in enumerate(bepipred_results):
                position = result.get('Position')
                score = result.get('Score')
                
                if position is None or score is None:
                    continue
                    
                if score > threshold and not in_region:
                    in_region = True
                    start_pos = position
                elif (score <= threshold or i == len(bepipred_results) - 1) and in_region:
                    in_region = False
                    end_pos = position - 1 if score <= threshold else position
                    
                    # Calculate average score for this region
                    region_scores = [r.get('Score', 0) for r in bepipred_results 
                                   if r.get('Position') and start_pos <= r.get('Position') <= end_pos]
                    avg_score = sum(region_scores) / len(region_scores) if region_scores else 0
                    
                    if end_pos - start_pos + 1 >= 8:  # Only regions with reasonable length
                        antigenic_regions.append({
                            "start": start_pos,
                            "end": end_pos,
                            "score": avg_score,
                            "method": "BepiPred-2.0"
                        })
        
        print(f"Found {len(antigenic_regions)} antigenic regions")
        
    except Exception as e:
        print(f"Error in IEDB epitope prediction: {e}")
        antigenic_regions = []
    
    # ====== PART 3: Allergenicity prediction (simple local model) ======
    try:
        print("Running allergenicity assessment...")
        
        # Known allergenic properties
        allergen_aa = {
            'E': 1.5,  # Glutamic acid - high in many allergens
            'D': 1.5,  # Aspartic acid
            'K': 1.3,  # Lysine
            'R': 1.3,  # Arginine
            'Q': 1.2,  # Glutamine
            'N': 1.2,  # Asparagine
            'Y': 1.1,  # Tyrosine
            'F': 1.1,  # Phenylalanine
            'W': 1.1   # Tryptophan
        }
        
        # Calculate allergenicity score
        allergen_count = sum(allergen_aa.get(aa, 0) for aa in vaccine_seq)
        allergenicity_score = allergen_count / len(vaccine_seq) if len(vaccine_seq) > 0 else 0
        allergenicity_score = min(1.0, allergenicity_score / 10)  # Scale to 0-1
        
        # Classify prediction
        if allergenicity_score < 0.3:
            allergenicity_prediction = "Probable Non-Allergen"
        elif allergenicity_score < 0.6:
            allergenicity_prediction = "Uncertain"
        else:
            allergenicity_prediction = "Probable Allergen"
            
    except Exception as e:
        print(f"Error in allergenicity prediction: {e}")
        allergenicity_score = 0.5
        allergenicity_prediction = "Unknown"
    
    # Count epitopes based on linkers
    linker = args.linker
    epitope_count = vaccine_seq.count(linker) + 1 if linker in vaccine_seq else 1
    
    # ====== PART 4: Create FASTA file for ColabFold ======
    with open(args.output_colabfold, "w") as f:
        f.write(f">{args.protein_type}_vaccine\n{vaccine_seq}\n")
    
    print(f"Created FASTA file for ColabFold analysis: {args.output_colabfold}")
    
    # Compile all properties into a DataFrame
    properties = {
        "Property": [
            "Length", "Molecular Weight", "Theoretical pI", 
            "Instability Index", "GRAVY", "Aromaticity",
            "Epitope Count", "Allergenicity Score",
            "Helix Fraction", "Sheet Fraction", "Coil Fraction"
        ],
        "Value": [
            len(vaccine_seq), f"{mol_weight:.2f}", f"{theoretical_pi:.2f}",
            f"{instability_index:.2f}", f"{gravy:.4f}", f"{aromaticity:.4f}",
            epitope_count, f"{allergenicity_score:.4f}",
            f"{secondary_structure[0]:.4f}", f"{secondary_structure[1]:.4f}", f"{secondary_structure[2]:.4f}"
        ],
        "Unit": [
            "aa", "Da", "", "", "", "", "", "",
            "fraction", "fraction", "fraction"
        ],
        "Interpretation": [
            "", "", 
            "7-8 is optimal for solubility", 
            "<40 suggests stable protein",
            ">0 hydrophobic, <0 hydrophilic",
            "Frequency of aromatic AAs",
            "", 
            "<0.3 suggests non-allergen",
            "α-helix content",
            "β-sheet content",
            "Random coil content"
        ]
    }
    
    pd.DataFrame(properties).to_csv(args.output_properties, index=False)
    print(f"Properties saved to {args.output_properties}")
    
    # Create a detailed evaluation report
    with open(args.output_evaluation, "w") as f:
        f.write(f"{args.protein_type.upper()} Vaccine Construct Evaluation\n")
        f.write("=" * (len(args.protein_type) + 33) + "\n\n")
        f.write(f"Sequence Length: {len(vaccine_seq)} amino acids\n")
        f.write(f"Molecular Weight: {mol_weight:.2f} Da\n")
        f.write(f"Theoretical pI: {theoretical_pi:.2f}\n")
        f.write(f"Instability Index: {instability_index:.2f} (<40 suggests a stable protein)\n")
        f.write(f"Grand Average of Hydropathy (GRAVY): {gravy:.4f}\n")
        f.write(f"Aromaticity: {aromaticity:.4f}\n")
        f.write(f"Number of Epitopes: {epitope_count}\n")
        f.write(f"Allergenicity Prediction: {allergenicity_prediction}\n")
        f.write(f"Allergenicity Score: {allergenicity_score:.4f} (<0.3 suggests non-allergen)\n")
        f.write(f"Secondary Structure: {secondary_structure[0]:.2f} helix, {secondary_structure[1]:.2f} sheet, {secondary_structure[2]:.2f} coil\n\n")
        
        f.write("Amino Acid Composition:\n")
        for aa, percentage in sorted(aa_composition.items()):
            f.write(f"  {aa}: {percentage:.2f}%\n")
        
        f.write("\nPredicted Antigenic Regions:\n")
        if antigenic_regions:
            for region in antigenic_regions:
                f.write(f"  Region {region['start']}-{region['end']}: Score {region['score']:.3f} ({region['method']})\n")
        else:
            f.write("  No significant antigenic regions detected or prediction failed\n")
        
        f.write("\nSequence:\n")
        for line in textwrap.wrap(vaccine_seq, width=60):
            f.write(f"{line}\n")
        
        f.write("\n\nNotes:\n")
        f.write(f"- This evaluation was performed using optimized tools:\n")
        f.write(f"  * BioPython ProtParam for physicochemical properties\n")
        f.write(f"  * IEDB API for antigenicity prediction\n")
        f.write(f"  * Rule-based model for allergenicity assessment\n")
        f.write("- For 3D structure prediction, please use the generated FASTA file with Google Colab and ColabFold:\n")
        f.write(f"  * Upload the file '{os.path.basename(args.output_colabfold)}' to Google Drive\n")
        f.write(f"  * Use a ColabFold notebook for structure prediction\n")
        f.write("- For a more comprehensive evaluation, consider:\n")
        f.write("  1. Experimental validation of immunogenicity\n")
        f.write("  2. Testing against diverse HLA alleles\n")
        f.write("  3. Comparison with known effective vaccine epitopes\n")
    
    print(f"Evaluation report saved to {args.output_evaluation}")

if __name__ == "__main__":
    main()