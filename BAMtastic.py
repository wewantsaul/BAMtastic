#!/usr/bin/env python3
"""
BAMtastic!
Runs samtools flagstat and mosdepth, generates CSV summary and visualization plots.
"""

import argparse
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
import sys
import re
import os

def run_samtools_flagstat(bam_file):
    """Run samtools flagstat and parse output."""
    try:
        result = subprocess.run(['samtools', 'flagstat', bam_file],
                              capture_output=True, text=True, check=True)
        return parse_flagstat_output(result.stdout, bam_file)
    except subprocess.CalledProcessError as e:
        print(f"Error running samtools flagstat on {bam_file}: {e}")
        return None
    except FileNotFoundError:
        print("Error: samtools not found. Please install samtools.")
        return None

def parse_flagstat_output(output, bam_file):
    """Parse samtools flagstat output."""
    lines = output.strip().split('\n')
    stats = {'sample': Path(bam_file).stem}

    for line in lines:
        if 'total' in line and 'secondary' not in line and 'supplementary' not in line:
            stats['total_reads'] = int(line.split()[0])
        elif 'mapped (' in line and 'primary' not in line:
            match = re.search(r'(\d+) \+ \d+ mapped \(([0-9.]+)%', line)
            if match:
                stats['mapped_reads'] = int(match.group(1))
                stats['mapping_rate'] = float(match.group(2))
        elif 'primary mapped' in line:
            match = re.search(r'(\d+) \+ \d+ primary mapped \(([0-9.]+)%', line)
            if match:
                stats['primary_mapped'] = int(match.group(1))
                stats['primary_mapping_rate'] = float(match.group(2))
        elif 'paired in sequencing' in line:
            stats['paired_reads'] = int(line.split()[0])
        elif 'properly paired' in line:
            match = re.search(r'(\d+) \+ \d+ properly paired \(([0-9.]+)%', line)
            if match:
                stats['properly_paired'] = int(match.group(1))
                stats['properly_paired_rate'] = float(match.group(2))
        elif 'duplicates' in line:
            stats['duplicates'] = int(line.split()[0])

    return stats

def run_mosdepth(bam_file, output_dir, window_size=500):
    """Run mosdepth for coverage analysis."""
    try:
        sample_name = Path(bam_file).stem
        prefix = Path(output_dir) / sample_name

        # Run mosdepth
        cmd = ['mosdepth', '-x', '-n', '--by', str(window_size), str(prefix), bam_file]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        return parse_mosdepth_output(str(prefix), sample_name)
    except subprocess.CalledProcessError as e:
        print(f"Error running mosdepth on {bam_file}: {e}")
        return None
    except FileNotFoundError:
        print("Error: mosdepth not found. Please install mosdepth.")
        return None

def parse_mosdepth_output(prefix, sample_name):
    """Parse mosdepth output files."""
    stats = {'sample': sample_name}

    # Parse summary file
    summary_file = f"{prefix}.mosdepth.summary.txt"
    if Path(summary_file).exists():
        df = pd.read_csv(summary_file, sep='\t')
        total_row = df[df['chrom'] == 'total']
        if not total_row.empty:
            stats['mean_depth'] = total_row['mean'].iloc[0]
            stats['min_depth'] = total_row['min'].iloc[0]
            stats['max_depth'] = total_row['max'].iloc[0]

    # Parse global distribution
    global_dist_file = f"{prefix}.mosdepth.global.dist.txt"
    if Path(global_dist_file).exists():
        df = pd.read_csv(global_dist_file, sep='\t', names=['depth', 'count'])
        total_bases = df['count'].sum()

        if total_bases == 0:
            print(f"Warning: mosdepth global distribution for sample '{sample_name}' has zero total bases. Skipping coverage stats.")
            stats['bases_with_coverage'] = 0
            stats['percent_covered'] = 0.0
            for threshold in [1, 5, 10, 20, 30]:
                stats[f'coverage_>=_{threshold}x'] = 0.0
        else:
            stats['bases_with_coverage'] = df[df['depth'] > 0]['count'].sum()
            stats['percent_covered'] = (stats['bases_with_coverage'] / total_bases) * 100

            # Coverage thresholds
            for threshold in [1, 5, 10, 20, 30]:
                covered = df[df['depth'] >= threshold]['count'].sum()
                stats[f'coverage_>=_{threshold}x'] = (covered / total_bases) * 100

    return stats

def combine_stats(flagstat_results, mosdepth_results):
    """Combine flagstat and mosdepth results."""
    combined = []

    # Create a dictionary for easy lookup
    mosdepth_dict = {result['sample']: result for result in mosdepth_results if result}

    for flagstat in flagstat_results:
        if flagstat is None:
            continue

        sample = flagstat['sample']
        combined_stats = flagstat.copy()

        # Add mosdepth stats if available
        if sample in mosdepth_dict:
            mosdepth_stats = mosdepth_dict[sample]
            combined_stats.update({k: v for k, v in mosdepth_stats.items() if k != 'sample'})

        combined.append(combined_stats)

    return combined

def create_summary_plots(df, output_dir):
    """Create summary plots."""
    plt.style.use('default')

    # Create figure with subplots
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('BAM Quality Control Summary', fontsize=16)

    # 1. Mapping Rate
    if 'mapping_rate' in df.columns:
        axes[0, 0].bar(df['sample'], df['mapping_rate'], color='skyblue')
        axes[0, 0].set_title('Mapping Rate (%)')
        axes[0, 0].set_ylabel('Percentage')
        axes[0, 0].tick_params(axis='x', rotation=45)

        # Add horizontal line at 80%
        axes[0, 0].axhline(y=80, color='red', linestyle='--', alpha=0.7, label='80% threshold')
        axes[0, 0].legend()

    # 2. Mean Depth
    if 'mean_depth' in df.columns:
        axes[0, 1].bar(df['sample'], df['mean_depth'], color='lightgreen')
        axes[0, 1].set_title('Mean Depth')
        axes[0, 1].set_ylabel('Depth')
        axes[0, 1].tick_params(axis='x', rotation=45)

    # 3. Percent Covered
    if 'percent_covered' in df.columns:
        axes[0, 2].bar(df['sample'], df['percent_covered'], color='coral')
        axes[0, 2].set_title('Genome Coverage (%)')
        axes[0, 2].set_ylabel('Percentage')
        axes[0, 2].tick_params(axis='x', rotation=45)

    # 4. Coverage Thresholds
    coverage_cols = [col for col in df.columns if 'coverage_>=_' in col]
    if coverage_cols:
        coverage_data = df[['sample'] + coverage_cols].set_index('sample')
        coverage_data.plot(kind='bar', ax=axes[1, 0], width=0.8)
        axes[1, 0].set_title('Coverage at Different Thresholds')
        axes[1, 0].set_ylabel('Percentage')
        axes[1, 0].tick_params(axis='x', rotation=45)
        axes[1, 0].legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    # 5. Total Reads
    if 'total_reads' in df.columns:
        axes[1, 1].bar(df['sample'], df['total_reads'], color='gold')
        axes[1, 1].set_title('Total Reads')
        axes[1, 1].set_ylabel('Number of Reads')
        axes[1, 1].tick_params(axis='x', rotation=45)

    # 6. Properly Paired Rate
    if 'properly_paired_rate' in df.columns:
        axes[1, 2].bar(df['sample'], df['properly_paired_rate'], color='mediumpurple')
        axes[1, 2].set_title('Properly Paired Rate (%)')
        axes[1, 2].set_ylabel('Percentage')
        axes[1, 2].tick_params(axis='x', rotation=45)

    plt.tight_layout()
    plt.savefig(Path(output_dir) / 'bam_qc_summary.png', dpi=300, bbox_inches='tight')
    plt.close()

def create_heatmap(df, output_dir):
    """Create heatmap of key metrics."""
    # Select numeric columns for heatmap
    numeric_cols = ['mapping_rate', 'mean_depth', 'percent_covered', 'properly_paired_rate']
    heatmap_cols = [col for col in numeric_cols if col in df.columns]

    if len(heatmap_cols) < 2:
        print("Not enough numeric columns for heatmap")
        return

    # Prepare data for heatmap
    heatmap_data = df[['sample'] + heatmap_cols].set_index('sample')

    # Normalize data for better visualization
    heatmap_normalized = heatmap_data.apply(lambda x: (x - x.min()) / (x.max() - x.min()) if x.max() > x.min() else x)

    plt.figure(figsize=(10, 8))
    sns.heatmap(heatmap_normalized.T, annot=True, fmt='.2f', cmap='RdYlBu_r',
                cbar_kws={'label': 'Normalized Values'})
    plt.title('BAM Quality Metrics Heatmap')
    plt.tight_layout()
    plt.savefig(Path(output_dir) / 'bam_qc_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser(
        description='Analyze BAM files using samtools flagstat and mosdepth',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  %(prog)s -i file1.bam file2.bam -o results/
  %(prog)s -i *.bam -o qc_results/ --window-size 1000
  %(prog)s -i sample.bam -o output/ --skip-mosdepth
        '''
    )

    parser.add_argument('-i', '--input', nargs='+', required=True,
                        help='Input BAM files')
    parser.add_argument('-o', '--output', required=True,
                        help='Output directory')
    parser.add_argument('--window-size', type=int, default=500,
                        help='Window size for mosdepth analysis (default: 500)')
    parser.add_argument('--skip-mosdepth', action='store_true',
                        help='Skip mosdepth analysis (only run flagstat)')
    parser.add_argument('--csv-only', action='store_true',
                        help='Only generate CSV, skip plots')

    args = parser.parse_args()

    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Validate input files
    bam_files = []
    for pattern in args.input:
        files = list(Path('.').glob(pattern))
        if not files:
            if Path(pattern).exists():
                bam_files.append(pattern)
            else:
                print(f"Warning: No files found matching {pattern}")
        else:
            bam_files.extend([str(f) for f in files])

    if not bam_files:
        print("Error: No valid BAM files found")
        sys.exit(1)

    print(f"Processing {len(bam_files)} BAM files...")

    # Run flagstat
    print("Running samtools flagstat...")
    flagstat_results = []
    for bam_file in bam_files:
        print(f"  Processing {bam_file}")
        result = run_samtools_flagstat(bam_file)
        flagstat_results.append(result)

    # Run mosdepth
    mosdepth_results = []
    if not args.skip_mosdepth:
        print("Running mosdepth...")
        for bam_file in bam_files:
            print(f"  Processing {bam_file}")
            result = run_mosdepth(bam_file, output_dir, args.window_size)
            mosdepth_results.append(result)

    # Combine results
    combined_stats = combine_stats(flagstat_results, mosdepth_results)

    if not combined_stats:
        print("Error: No valid results generated")
        sys.exit(1)

    # Create DataFrame
    df = pd.DataFrame(combined_stats)

    # Save CSV
    csv_file = output_dir / 'bam_qc_summary.csv'
    df.to_csv(csv_file, index=False)
    print(f"Summary saved to: {csv_file}")

    # Create plots
    if not args.csv_only:
        print("Creating plots...")
        create_summary_plots(df, output_dir)
        create_heatmap(df, output_dir)
        print(f"Plots saved to: {output_dir}")

    # Print summary
    print("\n" + "="*50)
    print("SUMMARY")
    print("="*50)

    if 'mapping_rate' in df.columns:
        print(f"Mapping rates: {df['mapping_rate'].min():.1f}% - {df['mapping_rate'].max():.1f}%")
    if 'mean_depth' in df.columns:
        print(f"Mean depths: {df['mean_depth'].min():.1f}x - {df['mean_depth'].max():.1f}x")
    if 'percent_covered' in df.columns:
        print(f"Genome coverage: {df['percent_covered'].min():.1f}% - {df['percent_covered'].max():.1f}%")

    # Flag problematic samples
    problems = []
    for _, row in df.iterrows():
        if 'mapping_rate' in row and row['mapping_rate'] < 50:
            problems.append(f"{row['sample']}: Low mapping rate ({row['mapping_rate']:.1f}%)")
        if 'percent_covered' in row and row['percent_covered'] < 50:
            problems.append(f"{row['sample']}: Low coverage ({row['percent_covered']:.1f}%)")

    if problems:
        print("\nPROBLEMATIC SAMPLES:")
        for problem in problems:
            print(f"  - {problem}")

    print(f"\nResults saved to: {output_dir}")

if __name__ == "__main__":
    main()
