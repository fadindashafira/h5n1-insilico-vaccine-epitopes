#!/usr/bin/env python3

import sys
import subprocess
import argparse
from datetime import datetime

def calculate_similarity(seq1, seq2):
    """Calculate sequence similarity between two peptides"""
    # For shorter sequences, check if one is contained in the other
    if len(seq1) <= len(seq2) and seq1 in seq2:
        return 1.0
    if len(seq2) <= len(seq1) and seq2 in seq1:
        return 1.0
    
    # For sequences of same length, calculate percent identity
    if len(seq1) == len(seq2):
        matches = sum(a == b for a, b in zip(seq1, seq2))
        return matches / len(seq1)
    
    # For different length sequences, use a sliding window approach
    min_len = min(len(seq1), len(seq2))
    max_similarity = 0
    
    # Check if shorter sequence is similar to any part of longer sequence
    if len(seq1) < len(seq2):
        shorter, longer = seq1, seq2
    else:
        shorter, longer = seq2, seq1
        
    for i in range(len(longer) - len(shorter) + 1):
        substring = longer[i:i+len(shorter)]
        matches = sum(a == b for a, b in zip(shorter, substring))
        similarity = matches / len(shorter)
        max_similarity = max(max_similarity, similarity)
    
    return max_similarity

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Design vaccine construct from epitopes')
    parser.add_argument('--combined-epitopes', required=True, help='Combined epitopes CSV file')
    parser.add_argument('--protein-type', required=True, help='Protein type (e.g., hemagglutinin, neuraminidase)')
    parser.add_argument('--linker', default='GPGPG', help='Linker sequence between epitopes')
    parser.add_argument('--max-epitopes', type=int, default=10, help='Maximum number of epitopes to include')
    parser.add_argument('--min-epitopes', type=int, default=1, help='Minimum number of epitopes required')
    parser.add_argument('--leading-seq', default='', help='Sequence to add at the N-terminus')
    parser.add_argument('--trailing-seq', default='', help='Sequence to add at the C-terminus')
    parser.add_argument('--output-fasta', required=True, help='Output FASTA file for vaccine construct')
    parser.add_argument('--output-report', required=True, help='Output HTML report file')
    parser.add_argument('--similarity-threshold', type=float, default=0.7, help='Similarity threshold for epitope diversity (0-1)')
    
    args = parser.parse_args()
    
    # Install required packages
    try:
        import pandas as pd
        import numpy as np
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        from Bio import SeqIO
    except ImportError:
        subprocess.check_call([sys.executable, '-m', 'pip', 'install', '--user', 'pandas', 'numpy', 'biopython'])
        import pandas as pd
        import numpy as np
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        from Bio import SeqIO
    
    # Define select_diverse_epitopes function inside main to access pandas
    def select_diverse_epitopes(epitopes_df, max_epitopes, similarity_threshold=0.7):
        """Select diverse set of epitopes based on score and sequence diversity"""
        selected = []
        
        # Sort by score descending
        sorted_epitopes = epitopes_df.sort_values('consensus_score', ascending=False)
        
        for _, epitope in sorted_epitopes.iterrows():
            seq = epitope['sequence']
            
            # Check if this sequence is too similar to already selected ones
            too_similar = False
            for selected_epitope in selected:
                similarity = calculate_similarity(seq, selected_epitope['sequence'])
                if similarity > similarity_threshold:
                    too_similar = True
                    break
                    
            # Only add if not too similar to existing selections
            if not too_similar:
                selected.append(epitope)
                
                # Stop once we have enough epitopes
                if len(selected) >= max_epitopes:
                    break
                    
        return pd.DataFrame(selected) if selected else pd.DataFrame(columns=epitopes_df.columns)
    
    import os
    
    print(f"Designing vaccine construct for {args.protein_type}")
    print(f"Using linker sequence: {args.linker}")
    print(f"Maximum number of epitopes: {args.max_epitopes}")
    print(f"Similarity threshold: {args.similarity_threshold}")
    
    # Read the combined epitopes
    try:
        epitopes_df = pd.read_csv(args.combined_epitopes)
        print(f"Loaded {len(epitopes_df)} epitopes from {args.combined_epitopes}")
    except Exception as e:
        print(f"Error reading epitope file: {e}")
        # Create empty dataframe if file is missing or empty
        epitopes_df = pd.DataFrame(columns=['sequence', 'type', 'consensus_score'])
    
    # Ensure we have some epitopes to work with
    if len(epitopes_df) == 0:
        print("Warning: No epitopes found in the input file.")
        epitope_sequences = []
        selected_epitopes = pd.DataFrame()
    else:
        # Ensure 'sequence' column exists
        if 'sequence' not in epitopes_df.columns:
            print("Error: 'sequence' column not found in epitopes file.")
            epitope_sequences = []
            selected_epitopes = pd.DataFrame()
        else:
            # Handle case where consensus_score might be missing
            if 'consensus_score' not in epitopes_df.columns:
                if 'score' in epitopes_df.columns:
                    epitopes_df['consensus_score'] = epitopes_df['score']
                else:
                    # Default score if missing
                    epitopes_df['consensus_score'] = 0.5
            
            # Ensure 'type' column exists
            if 'type' not in epitopes_df.columns:
                # Default to 'unknown' type
                epitopes_df['type'] = 'unknown'
            
            # Calculate target counts for each epitope type
            total_epitopes = min(args.max_epitopes, len(epitopes_df))
            
            # Two-step selection for diversity:
            # 1. First select diverse epitopes within each type
            b_cell_df = epitopes_df[epitopes_df['type'] == 'B-cell']
            mhc_i_df = epitopes_df[epitopes_df['type'] == 'MHC-I']
            mhc_ii_df = epitopes_df[epitopes_df['type'] == 'MHC-II']
            
            # Calculate proportional allocations
            b_cell_count = min(total_epitopes // 3, len(b_cell_df)) if len(b_cell_df) > 0 else 0
            mhc_i_count = min(total_epitopes // 3, len(mhc_i_df)) if len(mhc_i_df) > 0 else 0
            mhc_ii_count = min(total_epitopes - b_cell_count - mhc_i_count, len(mhc_ii_df)) if len(mhc_ii_df) > 0 else 0
            
            print(f"Allocating epitopes: B-cell: {b_cell_count}, MHC-I: {mhc_i_count}, MHC-II: {mhc_ii_count}")
            
            # Select diverse epitopes from each type
            bcell_epitopes = select_diverse_epitopes(b_cell_df, b_cell_count, args.similarity_threshold) if b_cell_count > 0 else pd.DataFrame(columns=epitopes_df.columns)
            mhci_epitopes = select_diverse_epitopes(mhc_i_df, mhc_i_count, args.similarity_threshold) if mhc_i_count > 0 else pd.DataFrame(columns=epitopes_df.columns)
            mhcii_epitopes = select_diverse_epitopes(mhc_ii_df, mhc_ii_count, args.similarity_threshold) if mhc_ii_count > 0 else pd.DataFrame(columns=epitopes_df.columns)
            
            # 2. Now combine and ensure diversity across the entire set
            combined_df = pd.concat([bcell_epitopes, mhci_epitopes, mhcii_epitopes])
            
            # If we don't have enough epitopes after the first round, try to fill from remaining
            if len(combined_df) < total_epitopes:
                remaining_count = total_epitopes - len(combined_df)
                print(f"Only selected {len(combined_df)} epitopes in first round, attempting to select {remaining_count} more")
                
                # Get sequences we've already selected
                selected_sequences = combined_df['sequence'].tolist()
                
                # Find epitopes we haven't selected yet
                remaining_df = epitopes_df[~epitopes_df['sequence'].isin(selected_sequences)]
                
                if len(remaining_df) > 0:
                    # Select additional diverse epitopes from remaining
                    additional_epitopes = select_diverse_epitopes(remaining_df, remaining_count, args.similarity_threshold)
                    combined_df = pd.concat([combined_df, additional_epitopes])
            
            # Final selection
            selected_epitopes = combined_df.sort_values('consensus_score', ascending=False).head(total_epitopes)
            
            # Extract sequences
            epitope_sequences = selected_epitopes['sequence'].tolist()
            
            print(f"Selected {len(epitope_sequences)} diverse epitopes for the vaccine construct")
    
    # Construct the vaccine sequence with linkers
    if len(epitope_sequences) >= args.min_epitopes:
        # Join epitopes with linkers
        middle_sequence = args.linker.join(epitope_sequences)
        
        # Add leading and trailing sequences if provided
        vaccine_sequence = args.leading_seq + middle_sequence + args.trailing_seq
    else:
        # If insufficient epitopes, log an error but still create a minimal valid output
        print(f"ERROR: Found only {len(epitope_sequences)} epitopes, which is less than the minimum required ({args.min_epitopes}).")
        
        # Create an empty but valid FASTA sequence - just the leading and trailing sequences if provided
        vaccine_sequence = args.leading_seq + args.trailing_seq
        
        if not vaccine_sequence:
            # If no leading/trailing sequences provided, use a single amino acid to create a valid FASTA
            vaccine_sequence = "M"
            print("WARNING: No epitopes, leading, or trailing sequences available. Creating minimal valid FASTA with just a Methionine.")
    
    # Create the vaccine record
    vaccine_record = SeqRecord(
        Seq(vaccine_sequence),
        id=f"{args.protein_type}_vaccine",
        name=f"{args.protein_type} epitope-based vaccine",
        description=f"H5N1 {args.protein_type} epitope-based vaccine construct designed on {datetime.now().strftime('%Y-%m-%d')}"
    )
    
    # Save the vaccine construct
    SeqIO.write(vaccine_record, args.output_fasta, "fasta")
    print(f"Vaccine construct saved to {args.output_fasta}")
    
    # Generate HTML report
    with open(args.output_report, 'w') as f:
        f.write(f'''<!DOCTYPE html>
    <html>
    <head>
        <title>{args.protein_type} Vaccine Design Report</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 20px; }}
            h1, h2 {{ color: #2c3e50; }}
            table {{ border-collapse: collapse; width: 100%; }}
            th, td {{ padding: 8px; text-align: left; border-bottom: 1px solid #ddd; }}
            th {{ background-color: #f2f2f2; }}
            .sequence {{ font-family: monospace; word-break: break-all; }}
            .stats {{ margin: 20px 0; padding: 10px; background-color: #f8f9fa; border-radius: 5px; }}
            .warning {{ color: #e74c3c; }}
        </style>
    </head>
    <body>
        <h1>H5N1 {args.protein_type.capitalize()} Epitope-Based Vaccine Design Report</h1>
        <p>Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        
        <div class="stats">
            <h2>Vaccine Statistics</h2>
            <p><strong>Protein:</strong> {args.protein_type}</p>
            <p><strong>Number of epitopes:</strong> {len(epitope_sequences)}</p>
            <p><strong>Linker sequence:</strong> {args.linker}</p>
            <p><strong>Total length:</strong> {len(vaccine_sequence)} amino acids</p>
            <p><strong>Similarity threshold:</strong> {args.similarity_threshold}</p>
            
            ''')
        
        # Add warning if insufficient epitopes
        if len(epitope_sequences) < args.min_epitopes:
            f.write(f'''
            <p class="warning"><strong>Warning:</strong> Only {len(epitope_sequences)} epitopes found, which is less than the minimum required ({args.min_epitopes}).</p>
            ''')
            
        f.write(f'''
        </div>
        
        <h2>Final Vaccine Construct</h2>
        <p class="sequence">{vaccine_sequence}</p>
        
        <h2>Construct Components</h2>
        <table>
            <tr>
                <th>Component</th>
                <th>Sequence</th>
                <th>Length</th>
            </tr>
        ''')
        
        # Add leading sequence if present
        if args.leading_seq:
            f.write(f'''
            <tr>
                <td>Leading sequence</td>
                <td class="sequence">{args.leading_seq}</td>
                <td>{len(args.leading_seq)}</td>
            </tr>
            ''')
            
        # Add epitopes section
        if epitope_sequences:
            f.write(f'''
            <tr>
                <td>Epitopes with linkers</td>
                <td class="sequence">{args.linker.join(epitope_sequences)}</td>
                <td>{len(args.linker.join(epitope_sequences))}</td>
            </tr>
            ''')
        
        # Add trailing sequence if present    
        if args.trailing_seq:
            f.write(f'''
            <tr>
                <td>Trailing sequence</td>
                <td class="sequence">{args.trailing_seq}</td>
                <td>{len(args.trailing_seq)}</td>
            </tr>
            ''')
            
        f.write('''
        </table>
        
        <h2>Selected Epitopes</h2>
        <table>
            <tr>
                <th>Sequence</th>
                <th>Type</th>
                <th>HLA Restriction</th>
                <th>Score</th>
            </tr>
        ''')
        
        # Add each epitope to the table
        if len(epitope_sequences) > 0:
            for _, row in selected_epitopes.iterrows():
                hla = row.get('hla', 'N/A')
                f.write(f'''
                <tr>
                    <td class="sequence">{row['sequence']}</td>
                    <td>{row['type']}</td>
                    <td>{hla}</td>
                    <td>{row.get('consensus_score', 'N/A'):.3f}</td>
                </tr>
                ''')
        else:
            f.write(f'''
            <tr>
                <td colspan="4">No epitopes selected.</td>
            </tr>
            ''')
        
        f.write('''
        </table>
    </body>
    </html>
        ''')
    
    print(f"Vaccine design report saved to {args.output_report}")

if __name__ == "__main__":
    main()