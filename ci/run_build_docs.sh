#!/bin/bash

echo "inside $0"

"$TRAVIS_BUILD_DIR"/ci/build_docs.sh 2>&1 > /tmp/doc.log &

# wait until subprocesses finish (build_docs.sh)
wait

exit 0
