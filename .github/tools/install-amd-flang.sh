#!/bin/sh

# Install amdflang
# 
#
# Originally written for Squash <https://github.com/quixdb/squash> by
# Evan Nemerson.  For documentation, bug reports, support requests,
# etc. please use <https://github.com/nemequ/pgi-travis>.
#
# To the extent possible under law, the author(s) of this script have
# waived all copyright and related or neighboring rights to this work.
# See <https://creativecommons.org/publicdomain/zero/1.0/> for
# details.

VERSION=7.0.5

TEMPORARY_FILES="${TMPDIR:-/tmp}"
export AMDFLANG_INSTALL_DIR=$(pwd)/amdflang-install
export AMDFLANG_SILENT=true
while [ $# != 0 ]; do
    case "$1" in
    "--prefix")
        export AMDFLANG_INSTALL_DIR="$2"; shift
        ;;
    "--tmpdir")
        TEMPORARY_FILES="$2"; shift
        ;;
    "--verbose")
        export AMDFLANG_SILENT=false;
        ;;
    "--version")
        VERSION="$2"; shift
        ;;
    *)
        echo "Unrecognized argument '$1'"
        exit 1
        ;;
    esac
    shift
done

# Example download URL for version 6.0.0
#    https://repo.radeon.com/rocm/misc/flang/rocm-afar-7450-drop-6.0.0-ubu.tar.bz2

ver="$(echo $VERSION | tr -d . )"
BASENAME=$(curl -s "https://repo.radeon.com/rocm/misc/flang/" | grep -oP "rocm-afar-[1-9][0-9]*-drop-{1}.{2}.{3}" | sort | tail -1)
URL_SHORT=$(curl -s "https://repo.radeon.com/rocm/misc/flang/" | grep -oP "$BASENAME-ubu[a-z]*.tar.bz2" | sort | tail -1)
URL=https://repo.radeon.com/rocm/misc/flang/${URL_SHORT}

if [ ! -f "${TEMPORARY_FILES}/${URL_SHORT}" ]; then
    echo "Downloading [${URL}]"
    wget -P ${TEMPORARY_FILES} "${URL}"
else
    echo "Download already present in ${TEMPORARY_FILES}"
fi

if [ ! -d "${AMDFLANG_INSTALL_DIR}/${BASENAME}" ]; then
    if [ ! -d "${AMDFLANG_INSTALL_DIR}" ]; then
        mkdir -p ${AMDFLANG_INSTALL_DIR}
    fi
    tar xjf ${TEMPORARY_FILES}/${URL_SHORT} -C ${AMDFLANG_INSTALL_DIR}
else
    echo "Install already present in ${AMDFLANG_INSTALL_DIR}"
fi

cat > ${AMDFLANG_INSTALL_DIR}/env.sh << EOF
### Variables
export AMDFLANG_INSTALL_DIR=${AMDFLANG_INSTALL_DIR}/${BASENAME}
export AMDFLANG_VERSION=${VERSION}

### Compilers
export PATH=\${AMDFLANG_INSTALL_DIR}/bin:\${PATH}
export LD_LIBRARY_PATH=\${AMDFLANG_INSTALL_DIR}/lib
export LD_LIBRARY_PATH=\${AMDFLANG_INSTALL_DIR}/llvm/lib:\$LD_LIBRARY_PATH
EOF

cat ${AMDFLANG_INSTALL_DIR}/env.sh