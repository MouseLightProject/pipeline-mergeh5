#! /usr/bin/env python

"""
Front-end script for merging an octree of .h5 files into a single big .h5 file.
"""

import argparse
import os
import subprocess
import sys

parser = argparse.ArgumentParser(description='Merge an octree of .h5 files into a single big .h5 file')
parser.add_argument('sample_date', help='Date of the sample to use, e.g. "2019-05-27"')
parser.add_argument('--threshold', nargs='?', help='Threshold to use, a 16-bit unsigned integer', default=6500)
args = parser.parse_args()

matlab_command_line_template = 'try; modpath; mergeh5_chunks(\'%s\', %d); catch err; fprintf(2, \'%%s\\n\', err.getReport()); quit(1); end; quit(0);'
print(matlab_command_line_template)
matlab_command_line = (matlab_command_line_template % (args.sample_date, args.threshold))
print(matlab_command_line)

script_file_path = os.path.abspath(__file__)
script_folder_path = os.path.dirname(script_file_path)
os.chdir(script_folder_path)
child = subprocess.Popen(['/misc/local/matlab-2018b/bin/matlab', '-nodisplay', '-r', matlab_command_line])
child.communicate()
rc = child.returncode
sys.exit(rc)
