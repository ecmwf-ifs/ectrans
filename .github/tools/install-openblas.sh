#! /usr/bin/env bash

# (C) Copyright 2025 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.

set +x
set -e -o pipefail

SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export PATH=$SCRIPTDIR:$PATH

os=$(uname)
case "$os" in
    Darwin)
      echo "Installing OpenBLAS via brew"
      brew ls --versions openblas || brew install openblas
      exit
    ;;
    Linux)
      echo "Installing OpenBLAS via apt-get"
      sudo apt-get install libblas-dev liblapack-dev
      exit
    ;;
    *)
    ;;
esac