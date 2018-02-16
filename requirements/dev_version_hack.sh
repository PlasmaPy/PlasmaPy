#!/bin/bash

# message for when hack starts executing
echo "================ Starting execution of dev_version_hack script ================="

# verbose printing of executed commands for debugging
if [[ $DEBUG == True ]]; then
    set -x
fi

# checking if PYTHON_VERSION variable is set
if [[ -z $PYTHON_VERSION ]]; then
    # fetches python --version output
    version="$(python --version 2>&1)"
    echo "Current version is $version"
    # grabbing just the number
    versionArr=($version)
    versionNum=${versionArr[1]}
    echo "Version num is $versionNum"
    # binds output of python --version to PYTHON_VERSION
    PYTHON_VERSION="$versionNum"
    echo "PYTHON_VERSION is now $PYTHON_VERSION"
fi

# disables debug if it was set
set +x

# message for when script is done
echo "================ Finished execution of dev_version_hack script ================="
