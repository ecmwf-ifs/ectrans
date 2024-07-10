#!/usr/bin/env bash

set -ea

ECTRANS_DOCS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/.."
cd ${ECTRANS_DOCS_DIR}

if [ -z "${ECTRANS_DOCS_TOKEN}" ]; then
    echo "Please provide credentials for https://sites.ecmwf.int/docs/ectrans/ "
    echo "  or provide a valid ECTRANS_DOCS_TOKEN in environment"
    echo -n "User: "
    read ECTRANS_DOCS_USER
    echo -n "Password: "
    read -s ECTRANS_DOCS_PASSWORD
    echo
fi

ford ford_config.md

# source venv/bin/activate

if [ -z "${ECTRANS_DOCS_TOKEN}" ]; then
    python scripts/publish.py --user=${ECTRANS_DOCS_USER} --password=${ECTRANS_DOCS_PASSWORD}
else
    python scripts/publish.py --token=${ECTRANS_DOCS_TOKEN}
fi
