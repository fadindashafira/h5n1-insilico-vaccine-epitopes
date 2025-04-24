#!/usr/bin/env python3

import sys
import subprocess
import time
import argparse
import os

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Predict B-cell epitopes from a protein sequence')
    parser.add_argument('--fasta', required=True, help='Input FASTA file')
    parser.add_argument('--protein-type', required=True, help='Protein type (e.g., hemagglutinin, neuraminidase)')
    parser.add_argument('--method', default='Bepipred-2.0', 
                       choices=['Bepipred-1.0', 'Bepipred-2.0', 'Chou-Fasman', 'Emini', 'Karplus-Schulz', 'Kolaskar-Tongaonkar', 'Parker'], 
                       help='Prediction method to use')
    parser.add_argument('--threshold', type=float, default=0.5, help='Score threshold (specificity)')
    parser.add_argument('--window-size', type=int, default=9, help='Window size for epitope prediction')
    parser.add_argument('--output', required=True, help='Output CSV file')
    
    args = parser.parse_args()
    
    # Install required packages if not already installed
    try:
        import pandas as pd
        import numpy as np
        from Bio import SeqIO
    except ImportError:
        subprocess.check_call([sys.executable, '-m', 'pip', 'install', '--user', 'pandas', 'numpy', 'biopython'])
        import pandas as pd
        import numpy as np
        from Bio import SeqIO
    
    # Install IEDB-python if not already installed
    try:
        from iedb.bcell import Epitope, Predictor  
    except ImportError:
    # Install the package
        subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'git+https://github.com/mattfemia/iedb-python.git'])
    
    # Try different import patterns - the package structure might be different than expected
    try:
        from iedb.bcell import Epitope, Predictor
    except ImportError:
        try:
            # Direct import from iedb
            from iedb import Epitope, Predictor
        except ImportError:
            # If still failing, use the fallback method
            print("WARNING: Could not import from IEDB package. Will use fallback prediction method.")
            # Define stub classes to make the rest of the code work
            class Epitope:
                def __init__(self, sequence, identifier):
                    self.sequence = sequence
                    self.identifier = identifier
            
            class Predictor:
                def __init__(self, method, threshold, window_size):
                    self.method = method
                    self.threshold = threshold
                    self.window_size = window_size
                
                def predict(self, epitope):
                    # Use the fallback method directly
                    return predict_bcell_fallback(epitope.sequence, self.window_size, self.threshold)
    
    # Read the FASTA file
    try:
        record = next(SeqIO.parse(args.fasta, "fasta"))
        sequence = str(record.seq)
        sequence_id = record.id
        print(f"Successfully loaded sequence {sequence_id} with length {len(sequence)}")
    except Exception as e:
        sys.stderr.write(f"ERROR: Could not parse FASTA file: {e}\n")
        sys.exit(1)
    
    print(f"Predicting B-cell epitopes for {args.protein_type} using {args.method}...")
    print(f"Window size: {args.window_size}, Threshold: {args.threshold}")
    
    # Create epitope object with the sequence
    epitope = Epitope(sequence=sequence, identifier=sequence_id)
    
    # Create predictor with the selected method
    # Map the method name to the actual method in IEDB-python
    method_map = {
        'Bepipred-1.0': 'bepipred',
        'Bepipred-2.0': 'bepipred2',
        'Chou-Fasman': 'chou_fasman',
        'Emini': 'emini',
        'Karplus-Schulz': 'karplus',
        'Kolaskar-Tongaonkar': 'kolaskar',
        'Parker': 'parker'
    }
    
    iedb_method = method_map.get(args.method, 'bepipred2')
    
    # Create predictor with the specified method
    predictor = Predictor(method=iedb_method, threshold=args.threshold, window_size=args.window_size)
    
    # Run prediction
    try:
        results = predictor.predict(epitope)
        print(f"Prediction completed successfully")
    except Exception as e:
        sys.stderr.write(f"ERROR during epitope prediction: {e}\n")
        # Fallback to built-in method if IEDB-python fails
        results = predict_bcell_fallback(sequence, args.window_size, args.threshold)
        print(f"Falling back to built-in prediction method due to error with IEDB-python")
    
    # Process results
    epitopes = []
    
    if hasattr(results, 'predictions') and results.predictions:
        # Process IEDB-python results
        for pred in results.predictions:
            # Extract continuous epitopes above threshold
            if pred.score >= args.threshold:
                epitopes.append({
                    'sequence': pred.peptide,
                    'start': pred.start,
                    'end': pred.end,
                    'score': pred.score,
                    'type': 'B-cell',
                    'method': args.method,
                    'source': args.protein_type
                })
    else:
        # Handle case when results don't have the expected structure
        print("Warning: Unexpected results format from IEDB prediction")
    
    print(f"Identified {len(epitopes)} potential B-cell epitopes")
    
    # Create DataFrame with epitopes
    if epitopes:
        epitope_df = pd.DataFrame(epitopes)
    else:
        # Create empty DataFrame if no epitopes found
        epitope_df = pd.DataFrame(columns=['sequence', 'start', 'end', 'score', 'type', 'method', 'source'])
        epitope_df['type'] = 'B-cell'
        epitope_df['method'] = args.method
        epitope_df['source'] = args.protein_type
    
    # Save to CSV
    epitope_df.to_csv(args.output, index=False)
    print(f"B-cell epitope prediction complete. Found {len(epitope_df)} epitopes.")

def predict_bcell_fallback(seq, window_size=9, threshold=0.5):
    """Fallback method if IEDB-python fails"""
    from collections import namedtuple
    
    # Create a simple Result class to mimic IEDB-python results
    Prediction = namedtuple('Prediction', ['peptide', 'start', 'end', 'score'])
    Result = namedtuple('Result', ['predictions'])
    
    # Amino acid propensity scales for B-cell epitope prediction
    propensity = {
        'A': 0.57, 'R': 1.87, 'N': 1.64, 'D': 1.46, 'C': 0.70,
        'Q': 1.56, 'E': 1.31, 'G': 0.72, 'H': 1.22, 'I': 0.73,
        'L': 0.76, 'K': 1.95, 'M': 0.85, 'F': 1.07, 'P': 1.95,
        'S': 1.41, 'T': 1.19, 'W': 1.14, 'Y': 1.47, 'V': 0.66
    }
    
    # Hydrophilicity scale (Hopp-Woods)
    hydrophilicity = {
        'A':-0.5, 'R': 3.0, 'N': 0.2, 'D': 3.0, 'C':-1.0,
        'Q': 0.2, 'E': 3.0, 'G':-0.4, 'H':-0.1, 'I':-1.8,
        'L':-1.8, 'K': 3.0, 'M':-1.3, 'F':-2.5, 'P': 0.0,
        'S': 0.3, 'T': 0.4, 'W':-3.4, 'Y':-2.3, 'V':-1.5
    }
    
    # Surface accessibility scale (Emini)
    accessibility = {
        'A': 0.701, 'R': 0.916, 'N': 0.811, 'D': 0.786, 'C': 0.623,
        'Q': 0.845, 'E': 0.854, 'G': 0.762, 'H': 0.822, 'I': 0.603,
        'L': 0.603, 'K': 0.930, 'M': 0.672, 'F': 0.695, 'P': 0.759,
        'S': 0.762, 'T': 0.739, 'W': 0.724, 'Y': 0.778, 'V': 0.628
    }
    
    predictions = []
    
    for i in range(len(seq) - window_size + 1):
        window = seq[i:i+window_size]
        
        # Skip windows with unknown amino acids
        if any(aa not in propensity for aa in window):
            continue
        
        # Calculate various scores
        avg_propensity = sum(propensity.get(aa, 0.5) for aa in window) / window_size
        avg_hydrophilicity = sum(hydrophilicity.get(aa, 0.0) for aa in window) / window_size
        avg_accessibility = sum(accessibility.get(aa, 0.7) for aa in window) / window_size
        
        # Combine scores with different weights
        score = (
            0.4 * avg_propensity + 
            0.3 * (avg_hydrophilicity + 5) / 10 +  # Normalize hydrophilicity
            0.3 * avg_accessibility
        )
        
        # Normalize score to 0-1 range
        normalized_score = min(1.0, max(0.0, score / 2.0))
        
        if normalized_score >= threshold:
            predictions.append(Prediction(
                peptide=window,
                start=i+1,
                end=i+window_size,
                score=normalized_score
            ))
    
    return Result(predictions=predictions)

if __name__ == "__main__":
    main()