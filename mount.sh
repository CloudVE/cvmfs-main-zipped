#!/bin/sh

mkdir -p /etc/cvmfs/
cp -r ./cvmfs-configs/ /etc/cvmfs/
mkdir -p /cvmfs/main/
mount -t cvmfs main.galaxyproject.org /cvmfs/main/


