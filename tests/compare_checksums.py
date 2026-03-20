#!/usr/bin/env python3

# (C) Copyright 2026- ECMWF.
# (C) Copyright 2026- RMI.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

import os
import filecmp
from glob import glob
import argparse

class colors:
    SUCCESS = '\033[94m'
    FAILURE = '\033[91m'
    ENDC = '\033[0m'

def compare_checksums(folder_path, ntasks, nthreads, exclude=""):
    if not os.path.isdir(folder_path):
        print(f"Error: '{folder_path}' is not a valid directory.")
        return False

    print(f"Check files in folder {folder_path}:")

    success_count = 0
    error_count = 0
    total_count = 0
    failed_list = []

    # Build list of all mpi0_omp1 checksum files
    reference_files = glob(os.path.join(folder_path, "*mpi0_omp1*.checksums"))
    if exclude:
        print(f"Excluding tests that match \"{exclude}\"")
        reference_files = [f for f in reference_files if exclude not in f]
    if len(reference_files) == 0:
        print("No reference checksum files found (mpi0_omp1)")
        return False

    # Search through all mpi0_omp1 reference checksum files
    for file_name in reference_files:
        if os.path.isfile(file_name):
            print(f"{file_name}")
            found = False
            for mpi in ntasks:
                for omp in nthreads:
                    other_file_name = file_name.replace("mpi0_omp1", f"mpi{mpi}_omp{omp}")
                    if other_file_name == file_name:
                        continue
                    if os.path.isfile(other_file_name):
                        total_count += 1
                        found = True
                        if (filecmp.cmp(file_name, other_file_name)):
                            print(f"    {other_file_name} ...{colors.SUCCESS}Passed{colors.ENDC}")
                            success_count += 1
                        else:
                            print(f"    {other_file_name} ...***{colors.FAILURE}Failed{colors.ENDC}")
                            error_count += 1
                            failed_list.append((file_name, other_file_name))
            if not found:
                print("No comparison found")
                return False
    percentage = int(100.0 * (success_count / total_count))
    if error_count > 0:
        print(f"{percentage}% comparisons passed, {colors.FAILURE}{error_count} comparisons failed out of "
              f"{total_count}{colors.ENDC}")

        print("The following checksum files do not match their mpi0_omp1 references:")

        for _, failed in failed_list:
            print(f"{colors.FAILURE}{failed}{colors.ENDC}")
    else:
        print(f"{percentage}% checks passed")

    return error_count == 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("folder_path", help="Path to the folder containing checksum files")
    parser.add_argument("ntasks", help="Comma-separated list of ntasks values to compare (e.g., '1,2,4')")
    parser.add_argument("nthreads", help="Comma-separated list of nthreads values to compare (e.g., '1,2,4')")
    parser.add_argument("-E", "--exclude", help="Pattern to exclude from comparison (e.g., gpu)", default="")
    args = parser.parse_args()

    folder = args.folder_path
    ntasks = args.ntasks.split(",")
    nthreads = args.nthreads.split(",")
    exclude = args.exclude

    if compare_checksums(folder, ntasks, nthreads, exclude=exclude):
        exit(0)
    else:
        exit(1)
