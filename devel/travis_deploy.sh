#!/bin/bash

set -e

# publish to pypi
# you must add encrypted TWINE_USERNAME and TWINE_PASSWORD env variables to .travis.yml
echo "//registry.npmjs.org/:_authToken=$NPM_TOKEN" > .npmrc
docker pull magland/mldevel_publish
docker run -e "TWINE_USERNAME=$TWINE_USERNAME" -e "TWINE_PASSWORD=$TWINE_PASSWORD" -v $PWD:/source magland/mldevel_publish pypi

