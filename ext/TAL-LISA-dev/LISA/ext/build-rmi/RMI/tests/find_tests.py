#! /usr/bin/env python3

import glob

for fn in glob.glob("*/Makefile"):
    print(fn.split("/")[0])
