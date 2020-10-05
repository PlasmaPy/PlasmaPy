#! /bin/bash

# Lists all the git commands required to push the changes to online repository.

git pull

# Remove any python compiled file.

find . -type f -name ".pyc" -exec rm -f {} \;

# Remove all the swap files. ( This might also remove recovery files. So 
# proceed with caution ).

find . -type f -name ".*.swp" -exec rm -f {} \;
find . -type f -name ".*.swo" -exec rm -f {} \;

# Add all the files which has been changed.

git add --a

# Commit all the added files.

git commit -a

# Push all the files to the GitHub

git push
