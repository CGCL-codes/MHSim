#!/bin/bash
export DEST="./.exvim.n"
export TOOLS="/home/xujh/exvim/vimfiles/tools/"
export TMP="${DEST}/_symbols"
export TARGET="${DEST}/symbols"
sh ${TOOLS}/shell/bash/update-symbols.sh
