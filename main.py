#!/usr/bin/env python

import sys
from mose.__main__ import run_with_sys_argv

argv = sys.argv[1:]
argv += ['--gui-mode', 'False', '--save-data', 'True', '--show-plot', 'True']

run_with_sys_argv(argv)
