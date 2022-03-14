export SOURCE_DIR=$HOME/source
export BUILD_DIR=$HOME/build
mkdir $SOURCE_DIR
mkdir $BUILD_DIR
cd $SOURCE_DIR

git clone https://github.com/sergeigribanov/ccgo
git clone https://github.com/sergeigribanov/KFBase
git clone https://github.com/sergeigribanov/KFCmd
git clone https://github.com/sergeigribanov/KFCmd5PiTest

source $PKG_DIR/root/bin/thisroot.sh

mkdir $BUILD_DIR/ccgo
cd $BUILD_DIR/ccgo
cmake -DCMAKE_INSTALL_PREFIX=$PKG_DIR/ccgo  -DCXX_STANDARD=17 $SOURCE_DIR/ccgo
make -j8
make install

mkdir $BUILD_DIR/KFBase
cd $BUILD_DIR/KFBase
cmake -DCMAKE_INSTALL_PREFIX=$PKG_DIR/KFBase -DCXX_STANDARD=17 $SOURCE_DIR/KFBase
make -j8
make install

mkdir $BUILD_DIR/KFCmd
cd $BUILD_DIR/KFCmd
cmake -DCMAKE_INSTALL_PREFIX=$PKG_DIR/KFCmd -DCXX_STANDARD=17 $SOURCE_DIR/KFCmd
make -j8
make install

mkdir $BUILD_DIR/KFCmd5PiTest
cd $BUILD_DIR/KFCmd5PiTest
cmake -DCMAKE_INSTALL_PREFIX=$PKG_DIR/KFCmd5PiTest -DCXX_STANDARD=17 $SOURCE_DIR/KFCmd5PiTest
make -j8
make install

cd $HOME/workdir
mv $HOME/yadisk.py .
python yadisk.py https://disk.yandex.ru/d/YKwqyx-GlnNRVQ -p $HOME/workdir
tar -xvf kfcmd_data.tar
