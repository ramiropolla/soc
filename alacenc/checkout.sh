echo "checking out ffmpeg svn"
rm -f ffmpeg/libavcodec/alacenc.c

svn checkout svn://svn.mplayerhq.hu/ffmpeg/trunk/ ffmpeg -r 13774
echo "patching ffmpeg"
cd ffmpeg
patch -p0 <../alacenc_build.patch
echo "copying the ALAC encoder source file to ffmpeg/libavcodec"
rm -f libavcodec/alacenc.c
ln -s ../../alacenc.c libavcodec/alacenc.c
echo "Done, now just do a regular configure and make to build."
