#!/bin/sh

git config --global user.email "travis@travis-ci.org"p
git config --global user.name "Travis CI"
git checkout master
git remote add origin-token https://almahmoud:${GH_TOKEN}@github.com/CloudVE/cvmfs-main-zipped.git
git add .
git commit -m "Travis update: $(date) (Build $TRAVIS_BUILD_NUMBER)" -m "[skip ci]" > ./current/commit.log
git commit -m "Commit Log of Travis update: $(date) (Build $TRAVIS_BUILD_NUMBER)" -m "[skip ci]"
git push --quiet origin-token master > /dev/null 2>&1
