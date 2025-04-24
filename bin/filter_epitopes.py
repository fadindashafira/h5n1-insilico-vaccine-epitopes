#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
from collections import defaultdict

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Filter and select epitopes for vaccine design')
    parser.add_argument('--bcell', required=True, help='B-cell epitopes CSV file')
    parser.add_argument('--tcell-i', required=True, help='T-cell MHC Class I epitopes CSV file')
    parser.add_argument('--tcell-ii', required=True, help='T-cell MHC Class II epitopes CSV file')
    parser.add_argument('--protein-type', required=True, help='Protein type (e.g., hemagglutinin, neuraminidase)')
    parser.add_argument('--output', required=True, help='Output CSV file for filtered epitopes')
    parser.add_argument('--min-score-bcell', type=float, default=0.6, help='Minimum score for B-cell epitopes')
    parser.add_argument('--min-score-tcell', type=float, default=0.6, help='Minimum score for T-cell epitopes')
    parser.add_argument('--overlap-threshold', type=int, default=5, help='Minimum overlap between epitopes for clustering')
    
    args = parser.parse_args()
    
    # Read epitope prediction results
    try:
        bcell_df = pd.read_csv(args.bcell)
        tcell_i_df = pd.read_csv(args.tcell_i)
        tcell_ii_df = pd.read_csv(args.tcell_ii)
        
        print(f"Loaded {len(bcell_df)} B-cell epitopes")
        print(f"Loaded {len(tcell_i_df)} T-cell MHC Class I epitopes")
        print(f"Loaded {len(tcell_ii_df)} T-cell MHC Class II epitopes")
    except Exception as e:
        print(f"Error loading epitope files: {e}")
        return
    
    # Filter epitopes by score
    bcell_filtered = bcell_df[bcell_df['score'] >= args.min_score_bcell]
    tcell_i_filtered = tcell_i_df[tcell_i_df['score'] >= args.min_score_tcell]
    tcell_ii_filtered = tcell_ii_df[tcell_ii_df['score'] >= args.min_score_tcell]
    
    print(f"After filtering by score:")
    print(f"  B-cell epitopes: {len(bcell_filtered)}")
    print(f"  T-cell MHC Class I epitopes: {len(tcell_i_filtered)}")
    print(f"  T-cell MHC Class II epitopes: {len(tcell_ii_filtered)}")
    
    # Function to check if two epitopes overlap
    def check_overlap(epitope1, epitope2, threshold=args.overlap_threshold):
        start1, end1 = epitope1['start'], epitope1['end']
        start2, end2 = epitope2['start'], epitope2['end']
        
        # Check if epitopes overlap
        if start1 <= end2 and start2 <= end1:
            overlap = min(end1, end2) - max(start1, start2) + 1
            return overlap >= threshold
        return False
    
    # Identify regions with both B-cell and T-cell epitopes
    combined_regions = []
    region_id = 1
    
    # Check each B-cell epitope for overlapping T-cell epitopes
    for _, b_epitope in bcell_filtered.iterrows():
        # Find overlapping MHC Class I epitopes
        overlapping_i = [t for _, t in tcell_i_filtered.iterrows() 
                         if check_overlap(b_epitope, t)]
        
        # Find overlapping MHC Class II epitopes
        overlapping_ii = [t for _, t in tcell_ii_filtered.iterrows() 
                          if check_overlap(b_epitope, t)]
        
        # If we have B-cell and at least one T-cell epitope overlapping
        if overlapping_i or overlapping_ii:
            # Determine the region span
            all_epitopes = [b_epitope] + overlapping_i + overlapping_ii
            start = min(e['start'] for e in all_epitopes)
            end = max(e['end'] for e in all_epitopes)
            
            # Calculate combined score (weighted average)
            b_score = b_epitope['score']
            ti_score = np.mean([t['score'] for t in overlapping_i]) if overlapping_i else 0
            tii_score = np.mean([t['score'] for t in overlapping_ii]) if overlapping_ii else 0
            
            # Weights for different epitope types
            w_b = 0.4  # B-cell weight
            w_ti = 0.3  # T-cell Class I weight
            w_tii = 0.3  # T-cell Class II weight
            
            # Calculate combined score with normalization
            if overlapping_i and overlapping_ii:
                combined_score = w_b * b_score + w_ti * ti_score + w_tii * tii_score
            elif overlapping_i:
                combined_score = (w_b * b_score + w_ti * ti_score) / (w_b + w_ti)
            elif overlapping_ii:
                combined_score = (w_b * b_score + w_tii * tii_score) / (w_b + w_tii)
            else:
                combined_score = b_score
            
            # Add to combined regions
            combined_regions.append({
                'region_id': region_id,
                'start': start,
                'end': end,
                'length': end - start + 1,
                'b_cell_score': b_score,
                'mhc_i_score': ti_score if overlapping_i else None,
                'mhc_ii_score': tii_score if overlapping_ii else None,
                'combined_score': combined_score,
                'b_cell_count': 1,
                'mhc_i_count': len(overlapping_i),
                'mhc_ii_count': len(overlapping_ii),
                'protein_type': args.protein_type
            })
            
            region_id += 1
    
    print(f"Identified {len(combined_regions)} regions with overlapping B and T cell epitopes")
    
    # Create DataFrame for combined regions
    if combined_regions:
        regions_df = pd.DataFrame(combined_regions)
        
        # Sort by combined score in descending order
        regions_df = regions_df.sort_values('combined_score', ascending=False)
        
        # Merge overlapping regions
        merged_regions = []
        processed = set()
        
        for i, region in regions_df.iterrows():
            if i in processed:
                continue
                
            current_region = region.to_dict()
            processed.add(i)
            
            # Check for overlapping regions to merge
            for j, other_region in regions_df.iterrows():
                if j in processed or i == j:
                    continue
                    
                # Check if regions overlap
                if (current_region['start'] <= other_region['end'] and 
                    other_region['start'] <= current_region['end']):
                    
                    # Merge regions
                    current_region['start'] = min(current_region['start'], other_region['start'])
                    current_region['end'] = max(current_region['end'], other_region['end'])
                    current_region['length'] = current_region['end'] - current_region['start'] + 1
                    
                    # Update counts and scores
                    current_region['b_cell_count'] += other_region['b_cell_count']
                    current_region['mhc_i_count'] += other_region['mhc_i_count']
                    current_region['mhc_ii_count'] += other_region['mhc_ii_count']
                    
                    # Take the max score as the merged score
                    current_region['combined_score'] = max(
                        current_region['combined_score'], 
                        other_region['combined_score']
                    )
                    
                    processed.add(j)
            
            merged_regions.append(current_region)
        
        print(f"After merging overlapping regions: {len(merged_regions)} regions")
        
        # Create final DataFrame with merged regions
        final_df = pd.DataFrame(merged_regions)
        
        # Sort by combined score in descending order
        final_df = final_df.sort_values('combined_score', ascending=False)
        
        # Save to CSV
        final_df.to_csv(args.output, index=False)
        print(f"Saved {len(final_df)} filtered epitope regions to {args.output}")
    else:
        # Create empty DataFrame if no regions found
        columns = ['region_id', 'start', 'end', 'length', 'b_cell_score', 
                   'mhc_i_score', 'mhc_ii_score', 'combined_score', 
                   'b_cell_count', 'mhc_i_count', 'mhc_ii_count', 'protein_type']
        
        final_df = pd.DataFrame(columns=columns)
        final_df.to_csv(args.output, index=False)
        print(f"No overlapping epitope regions found. Created empty output file.")

if __name__ == "__main__":
    main()