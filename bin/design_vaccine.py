#!/usr/bin/env python3

import sys
import subprocess
import argparse
from datetime import datetime

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
    
    import os
    
    print(f"Designing vaccine construct for {args.protein_type}")
    print(f"Using linker sequence: {args.linker}")
    print(f"Maximum number of epitopes: {args.max_epitopes}")
    
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
        # Group by type and select top epitopes from each
        b_cell_count = args.max_epitopes // 3
        mhc_i_count = args.max_epitopes // 3
        mhc_ii_count = args.max_epitopes - b_cell_count - mhc_i_count  # Ensure we use exactly max_epitopes
        
        bcell_epitopes = epitopes_df[epitopes_df['type'] == 'B-cell'].sort_values('consensus_score', ascending=False).head(b_cell_count)
        mhci_epitopes = epitopes_df[epitopes_df['type'] == 'MHC-I'].sort_values('consensus_score', ascending=False).head(mhc_i_count)
        mhcii_epitopes = epitopes_df[epitopes_df['type'] == 'MHC-II'].sort_values('consensus_score', ascending=False).head(mhc_ii_count)
        
        # Combine selected epitopes
        selected_epitopes = pd.concat([bcell_epitopes, mhci_epitopes, mhcii_epitopes])
        
        # Limit to max_epitopes
        selected_epitopes = selected_epitopes.sort_values('consensus_score', ascending=False).head(args.max_epitopes)
        
        # Extract sequences
        epitope_sequences = selected_epitopes['sequence'].tolist()
        
        print(f"Selected {len(epitope_sequences)} epitopes for the vaccine construct")
    
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