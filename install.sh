pip install -U requests

export SOURCE_DIR=$HOME/source
export BUILD_DIR=$HOME/build
mkdir $SOURCE_DIR
mkdir $BUILD_DIR
cd $SOURCE_DIR

git clone --branch develop https://github.com/sergeigribanov/KFBase
git clone --branch develop https://github.com/sergeigribanov/KFCmd
git clone https://github.com/sergeigribanov/gaussgen
git clone --branch develop https://github.com/sergeigribanov/KFCmd5PiTest

source $PKG_DIR/root/bin/thisroot.sh

mkdir $BUILD_DIR/KFBase
cd $BUILD_DIR/KFBase
cmake  -DCMAKE_CXX_STANDARD=17 -DCMAKE_INSTALL_PREFIX=$PKG_DIR/kfbase $SOURCE_DIR/KFBase
make -j8
make install

mkdir $BUILD_DIR/KFCmd
cd $BUILD_DIR/KFCmd
cmake  -DCMAKE_CXX_STANDARD=17 -DCMAKE_INSTALL_PREFIX=$PKG_DIR/kfcmd $SOURCE_DIR/KFCmd
make -j8
make install

mkdir $BUILD_DIR/gaussgen
cd $BUILD_DIR/gaussgen
cmake  -DCMAKE_CXX_STANDARD=17 -DCMAKE_INSTALL_PREFIX=$PKG_DIR/gaussgen $SOURCE_DIR/gaussgen
make -j8
make install

# mkdir $BUILD_DIR/KFCmd5PiTest
# cd $BUILD_DIR/KFCmd5PiTest
# cmake -DCMAKE_CXX_STANDARD=17 -DCMAKE_INSTALL_PREFIX=$PKG_DIR/KFCmd5PiTest $SOURCE_DIR/KFCmd5PiTest
# make -j8
# make install
