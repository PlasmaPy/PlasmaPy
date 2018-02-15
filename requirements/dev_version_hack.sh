#!/bin/bash
if [[ -z $PYTHON_VERSION ]]; then
    # fetches python --version output
    version="$(python --version 2>&1)"
    echo "Current version is $version"
    # binds output of python --version to PYTHON_VERSION
    PYTHON_VERSION="$version"
    echo "PYTHON_VERSION is now $PYTHON_VERSION"
fi
