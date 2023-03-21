pip install -U requests

export BUILD_DIR=$HOME/build
mkdir $BUILD_DIR
cd $SOURCE_DIR

git clone https://github.com/sergeigribanov/KFBase
git clone --branch tr_ph_v8 https://github.com/sergeigribanov/KFCmd KFCmd_tr_ph_v8
git clone --branch tr_ph_v9 https://github.com/sergeigribanov/KFCmd KFCmd_tr_ph_v9
git clone https://github.com/sergeigribanov/gaussgen

mkdir $BUILD_DIR/KFBase
cd $BUILD_DIR/KFBase
cmake  -DCMAKE_CXX_STANDARD=17 -DCMAKE_INSTALL_PREFIX=$PKG_DIR/kfbase $SOURCE_DIR/KFBase
make
make install

mkdir $BUILD_DIR/KFCmd_tr_ph_v8
cd $BUILD_DIR/KFCmd_tr_ph_v8
cmake  -DCMAKE_CXX_STANDARD=17 -DCMAKE_INSTALL_PREFIX=$PKG_DIR/kfcmd_tr_ph_v8 $SOURCE_DIR/KFCmd_tr_ph_v8
make
make install

mkdir $BUILD_DIR/KFCmd_tr_ph_v9
cd $BUILD_DIR/KFCmd_tr_ph_v9
cmake  -DCMAKE_CXX_STANDARD=17 -DCMAKE_INSTALL_PREFIX=$PKG_DIR/kfcmd_tr_ph_v9 $SOURCE_DIR/KFCmd_tr_ph_v9
make
make install

mkdir $BUILD_DIR/gaussgen
cd $BUILD_DIR/gaussgen
cmake  -DCMAKE_CXX_STANDARD=17 -DCMAKE_INSTALL_PREFIX=$PKG_DIR/gaussgen $SOURCE_DIR/gaussgen
make
make install
