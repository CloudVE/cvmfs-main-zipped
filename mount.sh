#!/bin/sh

mkdir -p /etc/
cp -r ./cvmfs/ /etc/
mkdir -p /cvmfs/main/
mount -t cvmfs main.galaxyproject.org /cvmfs/main/

