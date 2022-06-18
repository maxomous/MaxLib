#!/bin/bash


cd build;
if ! cmake .. ; then
    echo "CMake Failed"
	exit 1
fi

if ! make ; then
    echo "Make Failed"
	exit 1
fi


if ! sudo make install ; then
    echo "Install Failed"
	exit 1
fi

echo "Install Success"