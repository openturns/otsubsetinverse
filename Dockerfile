FROM ubuntu:16.04
RUN apt-get update && apt-get -y install cmake swig bison flex libmuparser-dev liblapack-dev libxml2-dev libboost-math-dev libtbb-dev python-dev python-scipy python-matplotlib python-numpydoc python-sphinx liblpsolve55-dev texlive-latex-recommended texlive-fonts-recommended texlive-latex-extra git
RUN git clone https://github.com/openturns/openturns.git && cd openturns && git checkout 1.9 && cmake -DCMAKE_INSTALL_PREFIX=$HOME/.local . && make install -j4 && cd ../
COPY . /usr/src/otsubsetinverse
WORKDIR /usr/src/otsubsetinverse/build
RUN rm -rf *
RUN cmake -DCMAKE_INSTALL_PREFIX=$HOME/.local -DUSE_SPHINX=OFF ..
RUN make -j10
RUN make install -j10