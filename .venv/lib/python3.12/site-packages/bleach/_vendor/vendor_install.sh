#!/bin/bash

set -e
set -u
set -o pipefail

BLEACH_VENDOR_DIR=${BLEACH_VENDOR_DIR:-"."}
DEST=${DEST:-"."}

# Install with no dependencies
pip install --no-binary all --no-compile --no-deps -r "${BLEACH_VENDOR_DIR}/vendor.txt" --target "${DEST}"

# Apply patches
(cd "${DEST}" && patch -p2 < 01_html5lib_six.patch)
(cd "${DEST}" && patch -p2 < 02_html5lib_wbr.patch)

# install Python 3.6.14 urllib.urlparse for #536
curl --proto '=https' --tlsv1.2 -o "${DEST}/parse.py" https://raw.githubusercontent.com/python/cpython/v3.6.14/Lib/urllib/parse.py
(cd "${DEST}" && sha256sum parse.py > parse.py.SHA256SUM)
