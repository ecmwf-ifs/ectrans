#!/usr/bin/env bash

set -ea

TOP_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../.." && pwd)"
INC_DIR="$TOP_DIR/src/trans/include/ectrans" # Location of ecTrans include files

# List include files in order to be processed
declare -a files=(
    "setup_trans0.h"
    "setup_trans.h"
    "trans_inq.h"
    "dir_trans.h"
    "inv_trans.h"
    "specnorm.h"
    "gpnorm_trans.h"
    "dist_grid.h"
    "gath_grid.h"
    "gath_spec.h"
    "get_current.h"
    "ini_spec_dist.h"
    "trans_pnm.h"
    "vordiv_to_uv.h"
    "dir_transad.h"
    "inv_transad.h"
    "trans_release.h"
    "trans_end.h"
    "dist_grid_32.h"
    "gath_grid_32.h"
)

# Prepare processed content directory
OUT_DIR=$TOP_DIR/docs/content_processed
rm -rf $OUT_DIR && mkdir $OUT_DIR
cp -r $TOP_DIR/docs/content/* $OUT_DIR

# Print list of all subroutines to API page
for file in "${files[@]}"; do
    # Get name of this subroutine in capitals
    subroutine_name=$(echo "${file%.*}" | tr '[:lower:]' '[:upper:]')

    # Print list of all subroutines
    echo "- [\`$subroutine_name\`](#$subroutine_name)" >> $OUT_DIR/api.md
    echo >> $OUT_DIR/api.md # Add a newline
done

# Build API documentation by extracting docblocks from every file in src/trans/include/ectrans
for file in "${files[@]}"; do
    # Get name of this subroutine in capitals
    subroutine_name=$(echo "${file%.*}" | tr '[:lower:]' '[:upper:]')

    # Extract docblock
    echo "<hr>" >> $OUT_DIR/api.md # Horizontal rule
    echo "<a name=\"$subroutine_name\"></a>" >> $OUT_DIR/api.md # Add an anchor
    awk '/begin_doc_block/{flag=1; next} /end_doc_block/{flag=0} flag' $INC_DIR/$file | \
        sed 's/^! //g' | sed 's/^!$//g' >> $OUT_DIR/api.md
    echo >> $OUT_DIR/api.md # Add a newline
done

# Generate FORD-based documentation of main ecTrans source code
ford ford_config.md

# Generate Doxygen-based documentation of transi
doxygen transi.doxygen
