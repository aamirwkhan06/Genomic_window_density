import argparse
import pandas as pd
import numpy as np
import os

def parse_coordinates(input_file):
    """
    Parse the input coordinate file into a pandas DataFrame
    """
    df = pd.read_csv(input_file, sep='\s+', header=None,
                     names=['chromosome', 'start', 'end', 'id', 'type'])
    return df

def parse_chr_lengths(chr_length_file):
    """
    Parse chromosome lengths file
    Expected format: chromosome<tab>length
    """
    chr_lengths = {}
    with open(chr_length_file, 'r') as f:
        for line in f:
            chrom, length = line.strip().split()
            chr_lengths[chrom] = int(length)
    return chr_lengths

def count_features_in_windows(coordinates_df, chr_lengths, window_size):
    """
    Count features in non-overlapping windows for each chromosome
    """
    results = []
    
    # Process each chromosome separately
    for chrom in sorted(chr_lengths.keys()):  # Sort chromosomes for consistent output
        chr_length = chr_lengths[chrom]
        
        # Get features for this chromosome
        chrom_data = coordinates_df[coordinates_df['chromosome'] == chrom]
        
        # Calculate number of full windows and remaining bases
        num_full_windows = chr_length // window_size
        remaining_bases = chr_length % window_size
        
        # Process full windows
        for i in range(num_full_windows):
            window_start = i * window_size
            window_end = (i + 1) * window_size
            
            # Count features that start in this window
            count = len(chrom_data[
                (chrom_data['start'] >= window_start) & 
                (chrom_data['start'] < window_end)
            ])
            
            results.append({
                'chromosome': chrom,
                'window_start': window_start,
                'window_end': window_end,
                'feature_count': count
            })
        
        # Process the last partial window if it exists
        if remaining_bases > 0:
            window_start = num_full_windows * window_size
            window_end = chr_length  # End at chromosome length
            
            count = len(chrom_data[
                (chrom_data['start'] >= window_start) & 
                (chrom_data['start'] < window_end)
            ])
            
            results.append({
                'chromosome': chrom,
                'window_start': window_start,
                'window_end': window_end,
                'feature_count': count
            })
    
    return pd.DataFrame(results)

def test_window_counting():
    """
    Test the window counting logic with a small dataset
    """
    print("Running tests...")
    
    # Create test data
    test_data = """Chr6\t1000\t2000\tID1\tTYPE1
Chr6\t2500\t3500\tID2\tTYPE1
Chr6\t9500\t11000\tID3\tTYPE1
Chr6\t9750\t10500\tID4\tTYPE1
Chr6\t15000\t16000\tID5\tTYPE1
Chr7\t500\t1500\tID6\tTYPE1
Chr7\t11000\t12000\tID7\tTYPE1"""
    
    # Write test data to temporary file
    with open('test_coords.txt', 'w') as f:
        f.write(test_data)
    
    # Create test chromosome lengths
    chr_lengths = {'Chr6': 20000, 'Chr7': 15000}
    with open('test_chr_lengths.txt', 'w') as f:
        for chrom, length in chr_lengths.items():
            f.write(f"{chrom}\t{length}\n")
    
    # Parse test data
    coords_df = parse_coordinates('test_coords.txt')
    
    # Test with 10kb windows
    window_size = 10000
    results = count_features_in_windows(coords_df, chr_lengths, window_size)
    
    # Print results for verification
    print("\nTest Results (10kb windows):")
    print(results)
    print("\nDetailed analysis of features in windows:")
    
    # Verify results
    for chrom, length in chr_lengths.items():
        chr_results = results[results['chromosome'] == chrom]
        
        # Check that the last window ends at chromosome length
        last_window = chr_results.iloc[-1]
        assert last_window['window_end'] == length, \
            f"Last window for {chrom} should end at {length}, but ends at {last_window['window_end']}"
        
        # Check that no window exceeds chromosome length
        assert all(results['window_end'] <= length), \
            f"Found windows extending beyond chromosome length for {chrom}"
    
    print("\nAll tests passed successfully!")
    
    # Cleanup test files
    os.remove('test_coords.txt')
    os.remove('test_chr_lengths.txt')

def main():
    parser = argparse.ArgumentParser(description='Count features in genomic windows')
    parser.add_argument('--test', action='store_true',
                       help='Run tests with sample data')
    parser.add_argument('input_file', nargs='?',
                       help='Input coordinates file')
    parser.add_argument('chr_lengths', nargs='?',
                       help='Chromosome lengths file')
    parser.add_argument('--window-size', type=int, default=100000,
                       help='Window size in base pairs (default: 100kb)')
    parser.add_argument('--output', default='window_counts.tsv',
                       help='Output file name (default: window_counts.tsv)')
    
    args = parser.parse_args()
    
    if args.test:
        test_window_counting()
        return
    
    if not args.input_file or not args.chr_lengths:
        parser.error("input_file and chr_lengths are required when not in test mode")
        
    # Parse input files
    coordinates_df = parse_coordinates(args.input_file)
    chr_lengths = parse_chr_lengths(args.chr_lengths)
    
    # Count features in windows
    results = count_features_in_windows(coordinates_df, chr_lengths, args.window_size)
    
    # Save results
    results.to_csv(args.output, sep='\t', index=False)
    print(f"Results saved to {args.output}")

if __name__ == "__main__":
    main()
