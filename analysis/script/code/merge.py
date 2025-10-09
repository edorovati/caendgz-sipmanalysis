import argparse
import os
import sys
import glob
# Add custom utils path relative to this file (../../../python)
base_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../python'))
sys.path.append(base_path)

from utils import Utils

# Argument parsing
parser = argparse.ArgumentParser(description="Merge and sum columns from multiple data files.")
parser.add_argument('--merge', type=str, nargs='+', required=True,
                    help='List of input files or a directory containing them (use wildcards or path).')
parser.add_argument('--columns', type=int, required=True, help='Total number of columns in each input file.')
parser.add_argument('--sum', dest='sum_columns', type=str, required=False, help='Columns to sum (e.g., 2,2-3,2-3,5)')
parser.add_argument('--output', type=str, required=True, help='Name of the output file.')
parser.add_argument('--mode', type=str, choices=['sum', 'append'], default='append',
                    help="Merge mode: 'sum' to add columns row-by-row, 'append' to concatenate rows.")

args = parser.parse_args()

# Parse and validate sum columns

if args.mode == 'sum':
    if not args.sum_columns:
        parser.error("--sum is required when --mode is 'sum'")
    sum_columns = Utils.parse_range(args.sum_columns)
    if max(sum_columns) >= args.columns:
        raise ValueError("Sum columns must be less than total columns specified.")
else:  # append mode
    if args.sum_columns:
        parser.error("--sum is not allowed when --mode is 'append'")
    sum_columns = []

# Run the merge operation
Utils.merge_files(args.merge, args.columns, sum_columns, args.output)
print(f"[INFO] File di output salvato in: {os.path.abspath(args.output)}")
