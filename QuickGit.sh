#!/bin/bash

if (( $# == 0 )) ; then
  	echo "Commit requires a Description"
	exit 1
fi	


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