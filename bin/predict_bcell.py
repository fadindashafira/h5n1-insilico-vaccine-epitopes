#!/usr/bin/env python3

import sys
import subprocess
import argparse
import os
import requests
import time
import pandas as pd
import numpy as np
from Bio import SeqIO

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Predict B-cell epitopes from a protein sequence')
    parser.add_argument('--fasta', required=True, help='Input FASTA file')
    parser.add_argument('--protein-type', required=True, help='Protein type (e.g., hemagglutinin, neuraminidase)')
    parser.add_argument('--method', default='Bepipred', 
                       choices=['Bepipred', 'Bepipred-2.0', 'Chou-Fasman', 'Emini', 'Karplus-Schulz', 'Kolaskar-Tongaonkar', 'Parker'], 
                       help='Prediction method to use')
    parser.add_argument('--threshold', type=float, default=0.5, help='Score threshold (specificity)')
    parser.add_argument('--window-size', type=int, default=9, help='Window size for epitope prediction')
    parser.add_argument('--output', required=True, help='Output CSV file')
    
    args = parser.parse_args()
    
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
    
    # Call IEDB API directly
    url = "http://tools-cluster-interface.iedb.org/tools_api/bcell/"
    
    # Prepare API request
    params = {
        'method': args.method,
        'sequence_text': sequence,
        'window_size': args.window_size
    }
    
    # Try API call with error handling
    max_retries = 3
    retry_delay = 5
    success = False
    
    for attempt in range(max_retries):
        try:
            print(f"Calling IEDB API (attempt {attempt+1}/{max_retries})...")
            response = requests.post(url, data=params, timeout=300)
            
            if response.status_code == 200:
                # Process successful response
                results = process_iedb_response(response.text, sequence, args.window_size, args.threshold)
                success = True
                break
            else:
                print(f"API call failed with status code {response.status_code}")
                print(f"Response: {response.text}")
        
        except Exception as e:
            print(f"Error during API call: {e}")
        
        # Wait before retry
        if attempt < max_retries - 1:
            print(f"Retrying in {retry_delay} seconds...")
            time.sleep(retry_delay)
    
    if not success:
        print("All API call attempts failed. Using fallback method...")
        results = predict_bcell_fallback(sequence, args.window_size, args.threshold)
    
    # Convert results to DataFrame
    if results:
        epitope_df = pd.DataFrame(results)
    else:
        # Create empty DataFrame if no epitopes found
        epitope_df = pd.DataFrame(columns=['sequence', 'start', 'end', 'score', 'type', 'method', 'source'])
    
    print(f"B-cell epitope prediction complete. Found {len(epitope_df)} epitopes.")
    
    # Save to CSV
    epitope_df.to_csv(args.output, index=False)

def process_iedb_response(response_text, sequence, window_size=9, threshold=0.5):
    """Process the IEDB API response to extract epitopes of specified length."""
    try:
        # Parse the response into a DataFrame for easier processing
        lines = response_text.strip().split('\n')
        
        # Skip the header line and parse the data
        data = []
        for line in lines[1:]:  # Skip header
            parts = line.split()
            if len(parts) >= 3:
                try:
                    position = int(parts[0])
                    residue = parts[1]
                    score = float(parts[2])
                    data.append({'position': position, 'residue': residue, 'score': score})
                except (ValueError, IndexError):
                    continue
        
        if not data:
            print("No data found in response")
            return []
        
        # Convert to DataFrame
        df = pd.DataFrame(data)
        
        # Calculate window averages to find epitopes of specified length
        epitopes = []
        for i in range(len(df) - window_size + 1):
            window_df = df.iloc[i:i+window_size]
            avg_score = window_df['score'].mean()
            
            if avg_score >= threshold:
                start_pos = int(window_df['position'].iloc[0])
                end_pos = int(window_df['position'].iloc[-1])
                # Extract sequence for this window
                epitope_seq = ''.join(window_df['residue'].tolist())
                
                epitopes.append({
                    'sequence': epitope_seq,
                    'start': start_pos,
                    'end': end_pos,
                    'score': avg_score,
                    'type': 'B-cell',
                    'method': 'Bepipred',
                    'source': 'IEDB API'
                })
        
        print(f"Identified {len(epitopes)} potential B-cell epitopes")
        return epitopes
    
    except Exception as e:
        print(f"Error processing API response: {e}")
        return []

def predict_bcell_fallback(seq, window_size=9, threshold=0.5):
    """Fallback method if IEDB API fails"""
    # Amino acid propensity scales for B-cell epitope prediction
    propensity = {
        'A': 0.57, 'R': 1.87, 'N': 1.64, 'D': 1.46, 'C': 0.70,
        'Q': 1.56, 'E': 1.31, 'G': 0.72, 'H': 1.22, 'I': 0.73,
        'L': 0.76, 'K': 1.95, 'M': 0.85, 'F': 1.07, 'P': 1.95,
        'S': 1.41, 'T': 1.19, 'W': 1.14, 'Y': 1.47, 'V': 0.66
    }
    
    results = []
    
    for i in range(len(seq) - window_size + 1):
        window = seq[i:i+window_size]
        
        # Skip windows with unknown amino acids
        if any(aa not in propensity for aa in window):
            continue
        
        # Calculate propensity score
        avg_propensity = sum(propensity.get(aa, 0.5) for aa in window) / window_size
        
        # Normalize score to 0-1 range
        normalized_score = min(1.0, max(0.0, avg_propensity / 2.0))
        
        if normalized_score >= threshold:
            results.append({
                'sequence': window,
                'start': i+1,
                'end': i+window_size,
                'score': normalized_score,
                'type': 'B-cell',
                'method': 'Fallback',
                'source': 'Local'
            })
    
    return results

if __name__ == "__main__":
    main()