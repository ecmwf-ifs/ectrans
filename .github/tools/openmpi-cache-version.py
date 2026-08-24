#!/usr/bin/env python3

import os
import subprocess
import sys


def main() -> int:
    try:
        prefix = subprocess.check_output(["brew", "--prefix", "open-mpi"], text=True).strip()
        resolved = os.path.realpath(prefix)
        if "/Cellar/" in resolved:
            print(os.path.basename(resolved))
            return 0
    except subprocess.CalledProcessError:
        pass

    try:
        version = subprocess.check_output(["mpirun", "--version"], text=True).splitlines()[0].split()[-1]
        print(version)
        return 0
    except subprocess.CalledProcessError:
        return 1


if __name__ == "__main__":
    sys.exit(main())