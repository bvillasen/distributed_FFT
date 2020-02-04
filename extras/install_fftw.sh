#!/bin/sh -e

myprefix=$HOME/apps/
FFTW_VERSION=3.3.5
INSTDIR=$myprefix/fftw-$FFTW_VERSION
TMP="${PWD}/tmp-fftw-$FFTW_VERSION"
LOGFILE="${TMP}/build.log"

bash check if directory exists
if [ -d $TMP ]; then
        echo "Directory $TMP already exists. Delete it? (y/n)"
	read answer
	if [ ${answer} = "y" ]; then
		rm -rf $TMP
	else
		echo "Program aborted."
		exit 1
	fi
fi

mkdir $TMP && cd $TMP
# cd $TMP
wget http://www.fftw.org/fftw-$FFTW_VERSION.tar.gz
tar -xzvf fftw-$FFTW_VERSION.tar.gz
cd fftw-$FFTW_VERSION

./configure --prefix=$INSTDIR -enable-mpi  \
  # CPPFLAGS=" -I/usr/lib/openmpi/include" \
  # LDFLAGS=" -L/usr/lib/openmpi/lib" \
  FC=f95 CC=gcc MPICC=mpiCC MPIFC=mpif90 2>&1 | tee $LOGFILE

make -j 4 2>&1 | tee -a $LOGFILE
make install 2>&1 | tee -a $LOGFILE
