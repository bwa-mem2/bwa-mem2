#! /usr/bin/env python3

import glob
import os

def all_tests():
    subfolders = [f.path for f in os.scandir(".") if f.is_dir()]
    tests = [x[2:] for x in subfolders if x != "./results"]
    return tests


any_failed = False
tests_remaining = set(all_tests())
for fn in sorted(glob.glob("results/*_result")):
    test_name = fn[8:-7]
    with open(fn) as f:
        first_line = int(next(f).strip())
        if first_line == 0:
            print("PASS", test_name)
        else:
            print("FAIL", test_name)
            any_failed = True

        tests_remaining.remove(fn[8:-7])

if tests_remaining:
    any_failed = True

    for missing in tests_remaining:
        print("MISS", missing)

if any_failed:
    exit(1)
exit(0)
