#!/bin/bash
export DEST="./.exvim.n"
export TOOLS="/home/xujh/exvim/vimfiles/tools/"
export TMP="${DEST}/_ID"
export TARGET="${DEST}/ID"
sh ${TOOLS}/shell/bash/update-idutils.sh
