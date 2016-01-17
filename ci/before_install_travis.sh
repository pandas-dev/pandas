#!/bin/bash

# If envars.sh determined we're running  in an authorized fork
# and the user opted in to the network cache,and that cached versions
# are available on the cache server, download and deploy the cached
# files to the local filesystem

echo "inside $0"

# overview
sudo apt-get update $APT_ARGS # run apt-get update for all versions

true # never fail because bad things happened here
