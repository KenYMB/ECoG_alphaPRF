#!/bin/bash
## Download data files from fileserver

# Get project directory
PROJECTNAME=ECoG_alpha
DATADIR=Data
if [ $# -ne 0 ];
then
    PROJECTDIR=$1 
else
    PROJECTDIR=$(cd "$(dirname $(readlink $0 || echo $0))/../..";pwd -P) 
fi
echo Dataset will be downloaded into ${PROJECTDIR}/${DATADIR}

# Set server directory
SERVERDIR=/Volumes/server/Projects/${PROJECTNAME}
if [ ! -e $SERVERDIR ]; then
    echo Type the path of winawerlab fileserver
    read SERVERDIR
    if [ -e ${SERVERDIR}/${PROJECTNAME} ]; then
        SERVERDIR=${SERVERDIR}/${PROJECTNAME}
    elif [ -e ${SERVERDIR}/${PROJECTNAME} ]; then
        SERVERDIR=${SERVERDIR}/Projects/${PROJECTNAME}
    elif [ -e ${SERVERDIR}/server/Projects/${PROJECTNAME} ]; then
        SERVERDIR=${SERVERDIR}/server/Projects/${PROJECTNAME}
    else
        echo "Couldn't find winawerlab fileserver."
        exit 1
    fi
fi

# Download
echo Copying from $SERVERDIR
rsync -avzuP ${SERVERDIR}/${DATADIR}  ${PROJECTDIR}
echo Complete
