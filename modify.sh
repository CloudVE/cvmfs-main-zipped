#!/bin/sh

rm -rf ./working/
mkdir ./working/
cp -r /cvmfs/main/config ./working/
cp -r /cvmfs/main/tool_data ./working/
find /cvmfs/main/shed_tools -type f -not -path "*/.hg/*" -path "*.xml" -exec cp --parents \{\} ./working/ \;
mv ./working/cvmfs/main/shed_tools ./working/shed_tools
rm -rf ./working/cvmfs/
find ./working/ -type f -exec sed -i.cvmfsoriginal -e 's%cvmfs/main.galaxyproject.org%cvmfs/zipped/cvmfs/main.galaxyproject.org%g' {} \;
rm -rf ./previous/
mv ./current/ ./previous/
mkdir -p ./current/original/
find ./working -name '*.cvmfsoriginal' -exec cp --parents \{\} ./current/ \;
mv ./current/working/ ./current/original
find ./working -name "*.cvmfsoriginal" -type f -delete
mkdir -p ./current/zipped/cvmfs
mv ./working/ main.galaxyproject.org
mv main.galaxyproject.org ./current/zipped/cvmfs
cd ./current
tar -zcvf ./zipped.tar.gz ./zipped/ > tar.log

