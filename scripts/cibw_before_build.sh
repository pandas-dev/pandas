#!/bin/bash
# Add 3rd party licenses, like numpy does
for file in $PACKAGE_DIR/LICENSES/*; do
  cat $file >> $PACKAGE_DIR/LICENSE
done
