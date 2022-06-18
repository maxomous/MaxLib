#!/bin/bash

# Use location of executable
cd $(dirname "$0")

if (( $# == 0 )) ; then
  	echo "Commit requires a Description"
	exit 1
fi	

# Add files (and clean up old)
if ! git add . ; then
    echo "Git Add Failed"
	exit 1
fi

if ! git commit -am "$1" ; then
    echo "Git Commit Failed"
	exit 1
fi


if ! git push ; then
    echo "Git Push Failed"
	exit 1
fi

echo "Git Success"