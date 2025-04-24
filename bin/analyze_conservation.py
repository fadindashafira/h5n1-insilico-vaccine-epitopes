#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
from Bio import AlignIO
from collections import Counter

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Analyze epitope conservation across H5N1 strains')
    parser.add_argument('--epitopes', required=True, help='Filtered epitopes CSV file')
    parser.add_argument('--alignment', required=True, help='Multiple sequence alignment file (FASTA, CLUSTAL, etc.)')
    parser.add_argument('--protein-type', required=True, help='Protein type (e.g., hemagglutinin, neuraminidase)')
    parser.add_argument('--threshold', type=float, default=0.9, help='Conservation threshold (0-1)')
    parser.add_argument('--output', required=True, help='Output CSV file for conserved epitopes')
    parser.add_argument('--format', default='fasta', help='Alignment file format')
    
    args = parser.parse_args()
    
    # Read epitope data
    try:
        epitopes_df = pd.read_csv(args.epitopes)
        print(f"Loaded {len(epitopes_df)} epitope regions")
    except Exception as e:
        print(f"Error loading epitope file: {e}")
        return
    
    # Read multiple sequence alignment
    try:
        alignment = AlignIO.read(args.alignment, args.format)
        print(f"Loaded alignment with {len(alignment)} sequences, each of length {alignment.get_alignment_length()}")
    except Exception as e:
        print(f"Error loading alignment file: {e}")
        return
    
    # Function to calculate conservation score for a position
    def calculate_position_conservation(position, alignment):
        # Extract all amino acids at this position
        column = [record.seq[position] for record in alignment]
        
        # Count occurrences of each amino acid
        counts = Counter(column)
        
        # Remove gaps from consideration
        if '-' in counts:
            gap_count = counts['-']
            del counts['-']
        else:
            gap_count = 0
        
        # Calculate conservation as frequency of most common amino acid
        if len(column) - gap_count > 0:
            most_common = counts.most_common(1)
            if most_common:
                return most_common[0][1] / (len(column) - gap_count)
            return 0
        return 0
    
    # Function to calculate conservation score for a region
    def calculate_region_conservation(start, end, alignment):
        scores = [calculate_position_conservation(pos, alignment) for pos in range(start-1, end)]
        return np.mean(scores) if scores else 0
    
    # Analyze conservation for each epitope region
    conserved_epitopes = []
    
    # Process epitopes with start/end positions
    if 'start' in epitopes_df.columns and 'end' in epitopes_df.columns:
        for _, epitope in epitopes_df.iterrows():
            # Make sure positions are within alignment range
            if epitope['start'] <= alignment.get_alignment_length() and epitope['end'] <= alignment.get_alignment_length():
                conservation_score = calculate_region_conservation(epitope['start'], epitope['end'], alignment)
                
                # Check if conservation is above threshold
                if conservation_score >= args.threshold:
                    # Create a new row with conservation data
                    epitope_dict = epitope.to_dict()
                    epitope_dict['conservation_score'] = conservation_score
                    conserved_epitopes.append(epitope_dict)
    
    # Process epitopes without positions but with sequences
    elif 'sequence' in epitopes_df.columns:
        # Function to find a sequence in alignment
        def find_sequence_positions(seq, reference_seq):
            seq = seq.replace('-', '')  # Remove gaps from query
            ref_str = str(reference_seq)
            
            # Find all occurrences of the sequence
            positions = []
            start = 0
            while True:
                start = ref_str.find(seq, start)
                if start == -1:
                    break
                positions.append((start + 1, start + len(seq)))
                start += 1
            
            return positions
            
        # For each epitope, try to find it in each alignment sequence
        for _, epitope in epitopes_df.iterrows():
            seq = epitope['sequence']
            conservation_counts = 0
            
            # Count how many sequences contain this epitope
            for record in alignment:
                positions = find_sequence_positions(seq, record.seq)
                if positions:
                    conservation_counts += 1
            
            # Calculate conservation rate
            conservation_score = conservation_counts / len(alignment) if len(alignment) > 0 else 0
            
            # Check if conservation is above threshold
            if conservation_score >= args.threshold:
                # Add to conserved epitopes
                epitope_dict = epitope.to_dict()
                epitope_dict['conservation_score'] = conservation_score
                conserved_epitopes.append(epitope_dict)
    
    print(f"Found {len(conserved_epitopes)} conserved epitopes with conservation >= {args.threshold}")
    
    # Create DataFrame and save to CSV
    if conserved_epitopes:
        conserved_df = pd.DataFrame(conserved_epitopes)
        conserved_df.to_csv(args.output, index=False)
        print(f"Saved conserved epitopes to {args.output}")
    else:
        # Create empty DataFrame with same columns as input plus conservation_score
        conserved_df = epitopes_df.copy()
        if len(conserved_df) > 0:
            conserved_df = conserved_df.head(0)  # Empty but keep columns
        conserved_df['conservation_score'] = None
        conserved_df.to_csv(args.output, index=False)
        print(f"No conserved epitopes found. Created empty output file.")

if __name__ == "__main__":
    main()