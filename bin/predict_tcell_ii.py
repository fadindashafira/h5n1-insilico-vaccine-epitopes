#!/usr/bin/env python3

import sys
import subprocess
import time
import argparse
import requests

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Predict T-cell MHC Class II epitopes from a protein sequence')
    parser.add_argument('--fasta', required=True, help='Input FASTA file')
    parser.add_argument('--protein-type', required=True, help='Protein type (e.g., hemagglutinin, neuraminidase)')
    parser.add_argument('--method', default='recommended', 
                       choices=['recommended', 'netmhciipan', 'nn_align', 'sturniolo', 'consensus3'],
                       help='Prediction method to use')
    parser.add_argument('--threshold', type=float, default=500, help='IC50 threshold (nM)')
    parser.add_argument('--length', type=int, default=15, help='Peptide length for MHC Class II epitopes')
    parser.add_argument('--alleles', help='Comma-separated list of HLA alleles')
    parser.add_argument('--output', required=True, help='Output CSV file')
    
    args = parser.parse_args()
    
    # Handle case insensitivity for method
    args.method = args.method.lower()
    
    # Install required packages
    try:
        import pandas as pd
        from Bio import SeqIO
    except ImportError:
        subprocess.check_call([sys.executable, '-m', 'pip', 'install', '--user', 'pandas', 'biopython', 'requests'])
        import pandas as pd
        from Bio import SeqIO
    
    # Read the FASTA file
    try:
        record = next(SeqIO.parse(args.fasta, "fasta"))
        sequence = str(record.seq)
        print(f"Successfully loaded sequence with length {len(sequence)}")
    except Exception as e:
        sys.stderr.write(f"ERROR: Could not parse FASTA file: {e}\n")
        sys.exit(1)
    
    # Parse alleles
    alleles = args.alleles.split(',') if args.alleles else []
    
    # If no alleles provided, use defaults
    if not alleles or not alleles[0]:
        alleles = ["HLA-DRB1*03:01", "HLA-DRB1*04:01", "HLA-DRB1*01:01", "HLA-DRB1*07:01"]
    
    print(f"Predicting T-cell MHC Class II epitopes for {args.protein_type} using {args.method}")
    print(f"Alleles: {alleles}")
    print(f"Peptide length: {args.length}")
    
    # IEDB MHC-II API URL
    url = "http://tools-cluster-interface.iedb.org/tools_api/mhcii/"
    
    # Process all alleles
    all_results = []
    
    for allele in alleles:
        print(f"Processing allele: {allele}")
        
        # Prepare API request data
        data = {
            'method': args.method,
            'sequence_text': sequence,
            'allele': allele,
            'length': args.length
        }
        
        # Make API call with retry mechanism
        max_retries = 3
        retry_delay = 5
        
        for attempt in range(max_retries):
            try:
                print(f"  Calling IEDB API (attempt {attempt+1}/{max_retries})...")
                response = requests.post(url, data=data, timeout=300)
                
                if response.status_code == 200:
                    try:
                        # Process the response
                        predictions = process_iedb_response(response.text, allele, args.threshold, sequence)
                        all_results.extend(predictions)
                        break
                    except Exception as e:
                        print(f"  Error processing response: {e}")
                else:
                    print(f"  API call failed with status code {response.status_code}")
                    print(f"  Response: {response.text[:500]}...")
            
            except Exception as e:
                print(f"  Error during API call: {e}")
            
            # Wait before retry
            if attempt < max_retries - 1:
                print(f"  Retrying in {retry_delay} seconds...")
                time.sleep(retry_delay)
    
    print(f"Total epitopes found: {len(all_results)}")
    
    # Create DataFrame with results
    if all_results:
        epitope_df = pd.DataFrame(all_results)
        epitope_df['type'] = 'MHC-II'
        epitope_df['method'] = args.method
        epitope_df['source'] = args.protein_type
        
        # Save to CSV
        epitope_df.to_csv(args.output, index=False)
        print(f"T-cell MHC-II epitope prediction complete. Found {len(epitope_df)} epitopes.")
    else:
        # Create empty DataFrame if no epitopes found
        epitope_df = pd.DataFrame(columns=['sequence', 'start', 'end', 'score', 'hla', 'ic50', 'percentile_rank', 'type', 'method', 'source'])
        epitope_df.to_csv(args.output, index=False)
        print(f"T-cell MHC-II epitope prediction complete. No epitopes found below threshold {args.threshold}.")

def process_iedb_response(response_text, allele, threshold, sequence):
    """Process the IEDB API response to extract epitopes."""
    results = []
    
    try:
        # Split into lines and skip header
        lines = response_text.strip().split('\n')
        if len(lines) <= 1:
            return results
            
        # Identify column positions from header
        header = lines[0].strip().split('\t')
        allele_col = header.index('allele') if 'allele' in header else 0
        peptide_col = header.index('peptide') if 'peptide' in header else 1
        ic50_col = -1
        rank_col = -1
        
        # Find IC50 and percentile rank columns
        for i, col in enumerate(header):
            if 'ic50' in col.lower():
                ic50_col = i
            elif 'rank' in col.lower():
                rank_col = i
        
        # Use default positions if columns not found
        if ic50_col == -1:
            # Class II output might have different column structure
            # Try to find any numeric column that might contain IC50 values
            for i in range(2, min(len(header), 6)):
                ic50_col = i
                break
                
        if rank_col == -1 and ic50_col != -1:
            rank_col = ic50_col + 1 if ic50_col + 1 < len(header) else ic50_col
        
        # Process each line (after header)
        for line in lines[1:]:
            parts = line.strip().split('\t')
            
            if len(parts) > max(peptide_col, ic50_col, 1):
                try:
                    # Extract data from response
                    peptide = parts[peptide_col] if peptide_col < len(parts) else "Unknown"
                    
                    # Get IC50 score - try to convert to float, use a high default if it fails
                    try:
                        ic50 = float(parts[ic50_col]) if ic50_col < len(parts) else float('inf')
                    except ValueError:
                        # Could be a percentile rank instead of IC50
                        ic50 = 5000
                    
                    # Get percentile rank if available
                    try:
                        percentile_rank = float(parts[rank_col]) if rank_col < len(parts) else 0
                    except ValueError:
                        percentile_rank = 0
                    
                    # Only include if below threshold
                    if ic50 <= threshold:
                        result = {
                            'sequence': peptide,
                            'start': 0,  # Will be updated if found in sequence
                            'end': 0,    # Will be updated if found in sequence
                            'score': 1.0 - (min(ic50, 5000) / 5000),  # Convert IC50 to 0-1 score
                            'hla': allele,
                            'ic50': ic50,
                            'percentile_rank': percentile_rank
                        }
                        
                        # Find position in sequence
                        pos = sequence.find(peptide)
                        if pos != -1:
                            result['start'] = pos + 1  # 1-based indexing
                            result['end'] = pos + len(peptide)
                            
                        results.append(result)
                except (ValueError, IndexError) as e:
                    print(f"  Error parsing line: {line}, error: {e}")
    
    except Exception as e:
        print(f"  Error processing response: {e}")
        
    return results

if __name__ == "__main__":
    main()