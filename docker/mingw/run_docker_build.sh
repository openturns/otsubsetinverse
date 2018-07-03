#!/bin/sh

set -x -e

uid=$1
gid=$2

cd /tmp

ARCH=x86_64
PYBASEVER=2.7
PYBASEVER_NODOT=${PYBASEVER:0:1}${PYBASEVER:2:2}
MINGW_PREFIX=/usr/${ARCH}-w64-mingw32/

export WINEARCH=win64
export WINEPATH=${MINGW_PREFIX}/bin
rm -rf ${HOME}/.wine*

export PYTHONHOME=${MINGW_PREFIX}
export PYTHONPATH=${MINGW_PREFIX}/lib/python${PYBASEVER_NODOT}
export MAKEFLAGS="-j8"

mkdir /tmp/build_otsubsetinverse && cd /tmp/build_otsubsetinverse
MOD_PREFIX=$PWD/install
MOD_PREFIX=${MINGW_PREFIX}
CXXFLAGS="-D_hypot=hypot -D_GLIBCXX_ASSERTIONS" ${ARCH}-w64-mingw32-cmake -DUSE_SPHINX=OFF \
  -DCMAKE_INSTALL_PREFIX=${MOD_PREFIX} \
  -DPYTHON_INCLUDE_DIR=/usr/${ARCH}-w64-mingw32/include/python${PYBASEVER_NODOT} \
  -DPYTHON_LIBRARY=/usr/${ARCH}-w64-mingw32/lib/libpython${PYBASEVER_NODOT}.dll.a \
  -DPYTHON_EXECUTABLE=/usr/${ARCH}-w64-mingw32/bin/python${PYBASEVER_NODOT}.exe \
  -DPYTHON_SITE_PACKAGES=Lib/site-packages \
  -DUSE_COTIRE=ON -DCOTIRE_MAXIMUM_NUMBER_OF_UNITY_INCLUDES="-j8" ../otsubsetinverse

make -j10
make install
make tests

# need some copies if install directory is not ${MINGW_PREFIX}
# cp ${MINGW_PREFIX}/bin/*.dll ${MOD_PREFIX}/bin
ctest  --output-on-failure --timeout 100 -j8

VERSION=`cat ../otsubsetinverse/VERSION`
cp /tmp/otsubsetinverse/distro/windows/* .
makensis -DMODULE_PREFIX=${MOD_PREFIX} -DMODULE_VERSION=${VERSION} -DOPENTURNS_VERSION=1.11 -DPYBASEVER=${PYBASEVER} -DPYBASEVER_NODOT=${PYBASEVER_NODOT} -DARCH=${ARCH} installer.nsi

# copy to host with same permission
if test -n "${uid}" -a -n "${gid}"
then
  sudo cp otsubsetinverse*.exe /tmp/otsubsetinverse
  sudo chown ${uid}:${gid} /tmp/otsubsetinverse/otsubsetinverse*.exe
fi
