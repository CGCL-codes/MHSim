#!/bin/bash
export DEST="./.exvim.n"
export TOOLS="/home/xujh/exvim/vimfiles/tools/"
export TMP="${DEST}/_inherits"
export TARGET="${DEST}/inherits"
sh ${TOOLS}/shell/bash/update-inherits.sh
