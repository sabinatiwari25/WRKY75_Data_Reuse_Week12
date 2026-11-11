"""Run the full analysis and render the submission after dependencies are installed."""
import argparse
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--rscript',default='Rscript',help='Path to Rscript, if it is not on PATH')
parser.add_argument('--fetch',action='store_true',help='Download and verify the GEO inputs first')
args = parser.parse_args()

def run(*command):
    subprocess.run(command,cwd=root,check=True)

if args.fetch:
    run(sys.executable,'scripts/fetch_geo.py')
run(sys.executable,'scripts/validate_inputs.py')
for name in ('check_software','analyze','verify_and_sensitivity'):
    run(args.rscript,f'scripts/{name}.R')
run(sys.executable,'scripts/render_report.py')
