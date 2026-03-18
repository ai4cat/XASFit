import argparse
from processor import XmuProcessor

'''Command-line interface for the xmu_extract.
   python xmu_extract.py \
    --src IN_PATH \
    --dst OUT_PATH \
    --marker "#  omega    e    k    mu    mu0     chi     @#"
'''

def parse_args():
    parser = argparse.ArgumentParser(
        description="Collect xmu.dat files and extract data after a marker."
    )
    parser.add_argument(
        '--src', type=str, required=True,
        help='Source folder containing xmu.dat files'
    )
    parser.add_argument(
        '--dst', type=str, required=True,
        help='Destination folder to save collected and extracted files'
    )
    parser.add_argument(
        '--marker', type=str, required=True,
        help='Marker line to locate data in .dat files'
    )
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    
    processor = XmuProcessor(args.src, args.dst, args.marker)
    print("copy done...")
    processor.collect()
    print("data extraction done...")
    processor.extract_data()
    print("Done!")
