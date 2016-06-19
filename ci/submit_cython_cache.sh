#!/bin/bash

CACHE_File="$HOME/.cache/cython_files.tar"
rm -rf $CACHE_File

home_dir=$(pwd)

pyx_files=`find ${TRAVIS_BUILD_DIR} -name "*.pyx"`
echo "pyx files:"
echo $pyx_files

tar cf ${CACHE_File} --files-from /dev/null

for i in ${pyx_files}
do
        f=${i%.pyx}
        ls $f.c* | tar rf  ${CACHE_File} -T -
done

echo "Cython files in cache tar:"
tar tvf ${CACHE_File}

exit 0
