echo "checking out ffmpeg svn"
svn checkout svn://svn.mplayerhq.hu/ffmpeg/trunk/ ffmpeg
echo "patching ffmpeg"
cd ffmpeg/libavcodec ;
patch -p0 <../../ffmpeg.patch
echo "copying the jpeg2000 files to ffmpeg/libavcodec"
ln -s ../../j2kenc.c j2kenc.c
ln -s ../../j2k.h j2k.h
echo "Done, now just do a regular configure and make to build."
