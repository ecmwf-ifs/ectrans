#!/usr/bin/env bash

set -ea

ECTRANS_DOCS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/.."
cd ${ECTRANS_DOCS_DIR}

if [ -z "${ECTRANS_DOCS_TOKEN}" ]; then
    echo "Please provide credentials for https://sites.ecmwf.int/docs/ectrans/ "
    echo "  by providing a valid ECTRANS_DOCS_TOKEN in environment"
    exit 1
fi

# Generate FORD-based documentation of main ecTrans source code
ford ford_config.md

# Generate Doxygen-based documentation of transi
doxygen transi.doxygen

# source venv/bin/activate

# Publish docs to sites.ecmwf.int/docs/ectrans
python scripts/publish.py --token=${ECTRANS_DOCS_TOKEN}
