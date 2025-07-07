#!/usr/bin/env python3

import sys
import os
import filecmp
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    
def compare_checksums(folder_path):
    if not os.path.isdir(folder_path):
        print(f"Error: '{folder_path}' is not a valid directory.")
        return
    
    print(f"Check files in folder {folder_path}:")
    
    success_count = 0
    error_count = 0
    total_count = 0
    failed_list = []
    for file_name in os.listdir(folder_path):        
        file_path = os.path.join(folder_path, file_name)
        if os.path.isfile(file_path) and "benchmark" in file_name and "mpi0_omp1" in file_name:
            print(f"{file_name}")
            found = False
            for mpi in [0,1,2]:
                for omp in [1,4,8]:                    
                    other_file_name = file_name.replace("mpi0_omp1",f"mpi{mpi}_omp{omp}")
                    if other_file_name == file_name:
                        continue
                    other_file_path = os.path.join(folder_path, other_file_name)
                    if os.path.isfile(other_file_path):
                        total_count = total_count + 1
                        found = True
                        if (filecmp.cmp(file_path, other_file_path)):
                            print(f"    {other_file_name} ...{bcolors.OKBLUE}Passed{bcolors.ENDC}")
                            success_count = success_count +1
                        else:
                            print(f"    {other_file_name} ...***{bcolors.FAIL}Failed{bcolors.ENDC}")
                            error_count = error_count + 1
                            failed_list.append(f"{file_path} {other_file_path}")
            if (not found):
                print(f"    No comparison found")
    percentage = int(100*(success_count/total_count))
    if (error_count> 0):
        print(f"{percentage}% comparison passed, {bcolors.FAIL}{error_count} comparison failed out of {total_count}{bcolors.ENDC}")
    
        print("The following comparisons FAILED:")
        
        for failed in failed_list:
            print(f"    {bcolors.FAIL}{failed}{bcolors.ENDC}")
    else:
        print(f"{percentage}% checks passed")

    if (error_count > 0):
        return False
    return True

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python compare_checksums.py <folder_path>")
    else:
        folder = sys.argv[1]
        result = compare_checksums(folder)
