"""
Main module. Run 'python -m mose' at the directory that contains the
directory 'mose'.
"""
import sys
from run import run_with_sys_argv

# Set options when running from a command prompt
run_with_sys_argv(sys.argv[1:])
