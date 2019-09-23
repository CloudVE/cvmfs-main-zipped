#!/bin/sh

sudo rm -rf ./working/
sudo mkdir ./working/
sudo cp -r /cvmfs/main/config ./working/
sudo cp -r /cvmfs/main/tool_data ./working/
sudo find ./working/ -type f -exec sed -i.cvmfsoriginal -e 's%cvmfs/main.galaxyproject.org%cvmfs/zipped%g' {} \;
sudo rm -rf ./previous/
sudo mv ./current/ ./previous/
sudo mkdir -p ./current/original/
find ./working -name '*.cvmfsoriginal' -exec cp --parents \{\} ./current/original \;
find ./working -name "*.cvmfsoriginal" -type f -delete
sudo mv ./working/ ./current/zipped/
sudo tar -zcvf ./current/zipped.tar.gz ./current/zipped/

