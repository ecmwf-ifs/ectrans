#!/usr/bin/env python3

import sys
import os
import filecmp
class colors:
    SUCCESS = '\033[94m'
    FAILURE = '\033[91m'
    ENDC = '\033[0m'

def compare_checksums(folder_path, ntasks, nthreads):
    if not os.path.isdir(folder_path):
        print(f"Error: '{folder_path}' is not a valid directory.")
        return False

    print(f"Check files in folder {folder_path}:")

    success_count = 0
    error_count = 0
    total_count = 0
    failed_list = []
    for file_name in os.listdir(folder_path):
        if (".checksums" in file_name and "benchmark" in file_name  and "mpi0_omp1" in file_name):
            file_path = os.path.join(folder_path, file_name)
            if os.path.isfile(file_path):
                print(f"{file_name}")
                found = False
                for mpi in ntasks:
                    for omp in nthreads:
                        other_file_name = file_name.replace("mpi0_omp1",f"mpi{mpi}_omp{omp}")
                        if other_file_name == file_name:
                            continue
                        other_file_path = os.path.join(folder_path, other_file_name)
                        if os.path.isfile(other_file_path):
                            total_count = total_count + 1
                            found = True
                            if (filecmp.cmp(file_path, other_file_path)):
                                print(f"    {other_file_name} ...{colors.SUCCESS}Passed{colors.ENDC}")
                                success_count = success_count +1
                            else:
                                print(f"    {other_file_name} ...***{colors.FAILURE}Failed{colors.ENDC}")
                                error_count = error_count + 1
                                failed_list.append(f"{file_path} {other_file_path}")
                if (not found):
                    print(f"    No comparison found")
                    return False
    percentage = int(100*(success_count/total_count))
    if (error_count> 0):
        print(f"{percentage}% comparison passed, {colors.FAILURE}{error_count} comparison failed out of {total_count}{colors.ENDC}")

        print("The following comparisons FAILED:")

        for failed in failed_list:
            print(f"    {colors.FAILURE}{failed}{colors.ENDC}")
    else:
        print(f"{percentage}% checks passed")

    if (error_count > 0):
        return False
    return True

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python compare_checksums.py <folder_path> <ntasks list> <nthreads list>")
        exit(1)
    else:
        folder = sys.argv[1]
        ntasks = sys.argv[2].split(",")
        nthreads = sys.argv[3].split(",")
        if compare_checksums(folder, ntasks, nthreads):
            exit(0)
        else:
            exit(1)
