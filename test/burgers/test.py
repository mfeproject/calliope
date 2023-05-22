#!/usr/bin/env python3

import numpy as np
import sys

def compare_max(field1, field2, tol, name, frame):
    """Calculates the maximum absolute difference between two fields,
    and compares against an input tolerance. Returns 0 if error is below
    tolerance and 1 otherwise. Prints result of comparison to terminal.
    Fields may be a scalar value."""
    err = abs(field1 - field2).max()
    report_test(frame, name, err, tol, "max")
    return 1 if not err <= tol else 0 # checks if err is nan

def report_test(frame, name, err, tol, compare_type):
    status = "PASS" if err <= tol else "FAIL"
    print("{:s}: frame {:d} {:s} {:s} error={:8.2e} (tol={:8.2e})"
          .format(status, frame, name, compare_type, err, tol))


nfail = 0
rows_per_frame = 42

for f in [1, 2, 5, 15, 28]:
  skiprows = 1 + rows_per_frame*(f-1)

  data = np.loadtxt("mfegrf", skiprows=skiprows, max_rows=rows_per_frame-1)
  gold = np.loadtxt(sys.argv[1], skiprows=skiprows, max_rows=rows_per_frame-1)

  nfail += compare_max(data[:,0], gold[:,0], tol=1e-7, name="x", frame=f)
  nfail += compare_max(data[:,1], gold[:,1], tol=1e-7, name="u", frame=f)

assert nfail == 0
