cd cvodes
tar -xzf cvodes-2.7.0.tar.gz
cd cvodes-2.7.0
./configure --prefix=$PWD/cvodes
make
make install
cd ../
cd ../
tar -xzf StochKit2.0.10.tgz
cd StochKit2.0.10
./install.sh
cd ../
export STOCHKIT_HOME=$PWD/StochKit2.0.10
export STOCHKIT_ODE=$PWD
make
