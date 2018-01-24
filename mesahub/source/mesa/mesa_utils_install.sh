set -x 
cd /tmp/mesa
tar -zxvf ndiff-2.00.tar.gz 
cd ndiff-2.00
./configure 
make all
make install
cd ..
tar -xzvf makedepf90-2.8.8.tar.gz 
cd makedepf90-2.8.8
./configure
make
make install

