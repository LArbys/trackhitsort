#!/bin/bash

home=$PWD

cd ../dllee_unified/
source configure.sh

cd $home
export LD_LIBRARY_PATH=$home:${LD_LIBRARY_PATH}
