#!/bin/sh

sudo mkdir -p /etc/cvmfs/
sudo cp -r ./cvmfs-configs/ /etc/cvmfs/
sudo mkdir -p /cvmfs/main/
sudo mount -t cvmfs main.galaxyproject.org /cvmfs/main/


