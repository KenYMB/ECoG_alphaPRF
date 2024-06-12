#!/bin/bash
## Download data files from OpenNeuro
#   ecog_APRFF_00_downloadData.sh DIRECTORY

## Check openneuro function
if ! type openneuro >/dev/null 2>&1; then
    echo "(Error) Need to install openneuro in your system" >&2
    exit 1
fi

# Get project directory
if [ $# -ne 0 ];
then
    DATADIR=$1 
else
    DATADIR=$(cd "$(dirname $(readlink $0 || echo $0))/..";pwd -P)/BIDS 
fi
echo Dataset will be downloaded into ${DATADIR}
mkdir -p ${DATADIR}

# OpenNeuro
ACCESSID=ds004194
VERSION=1.0.1
openneuro download -s $VERSION $ACCESSID ${DATADIR}

# Download
echo Complete
