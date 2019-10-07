#!/bin/sh

git config --global user.email "travis@travis-ci.org"p
git config --global user.name "Travis CI"
git add .
git add ./current
git add ./previous
git commit -m "Travis update: $(date) (Build $TRAVIS_BUILD_NUMBER)" -m "[skip ci]" > ./current/commit.log
git add ./current
git commit -m "Commit Log of Travis update: $(date) (Build $TRAVIS_BUILD_NUMBER)" -m "[skip ci]"
git remote add origin-token https://almahmoud:${GH_TOKEN}@github.com/CloudVE/cvmfs-main-zipped.git
git push --quiet origin-token master > /dev/null 2>&1
