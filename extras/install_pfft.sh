#!/bin/sh -e

myprefix=$HOME/apps/
PFFT_VERSION=git
FFTW_VERSION=3.3.5
INSTDIR=$myprefix/pfft-$PFFT_VERSION
FFTWDIR=/home/bruno/apps/fftw-3.3.5
TMP="${PWD}/tmp-pfft-$PFFT_VERSION"
LOGFILE="${TMP}/build.log"

# #TORNADO
# MPIDIR=/usr/lib/openmpi

#HAWK
MPIDIR=/usr/lib/x86_64-linux-gnu/openmpi

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
#wget https://github.com/mpip/pfft/archive/master.tar.gz
git clone https://github.com/mpip/pfft.git
# wget http://www.tu-chemnitz.de/~potts/workgroup/pippig/software/pfft-$PFFT_VERSION.tar.gz
# gzip -dc master | tar xvf -
# cd pfft-$PFFT_VERSION
# cp ./tmp-pfft-1.0.8-alpha/pfft-1.0.8-alpha/configure ./tmp-pfft-git/pfft-master/
cd pfft
./bootstrap.sh
# ./configure --prefix=$INSTDIR --disable-shared --enable-openmp \
./configure --prefix=$INSTDIR --disable-shared  \
  CPPFLAGS="-I$FFTWDIR/include -I$MPIDIR/include" \
  LDFLAGS="-L$FFTWDIR/lib -L$MPIDIR/lib" \
  FC=f95 CC=gcc MPICC=mpiCC MPIFC=mpif90 2>&1 | tee $LOGFILE

make -j 4 2>&1 | tee -a $LOGFILE
make install 2>&1 | tee -a $LOGFILE
