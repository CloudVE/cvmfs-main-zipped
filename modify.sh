#!/bin/sh

rm -rf ./working/
mkdir ./working/
cp -r /cvmfs/main/config ./working/
cp -r /cvmfs/main/tool_data ./working/
find /cvmfs/main/shed_tools -type f -not -path "*/test-data/*" -not -path "*/.hg/*" -exec cp --parents \{\} ./working/ \;
mv ./working/cvmfs/main/shed_tools ./working/shed_tools
rm -rf ./working/cvmfs/main/
find ./working/ -type f -exec sed -i.cvmfsoriginal -e 's%cvmfs/main.galaxyproject.org%cvmfs/zipped%g' {} \;
rm -rf ./previous/
mv ./current/ ./previous/
mkdir -p ./current/original/
find ./working -name '*.cvmfsoriginal' -exec cp --parents \{\} ./current/ \;
mv ./current/working/ ./current/original
find ./working -name "*.cvmfsoriginal" -type f -delete
mv ./working/ ./current/zipped/
tar -zcvf ./current/zipped.tar.gz ./current/zipped/ > ./current/tar.log

