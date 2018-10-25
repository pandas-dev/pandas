#!/usr/bin/env bash

echo "[installed versions]"

export PATH="$MINICONDA_DIR/bin:$PATH"
source activate pandas

python -c "import pandas; pandas.show_versions();"
