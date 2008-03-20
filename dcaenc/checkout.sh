echo "checking out ffmpeg svn"
svn checkout svn://svn.mplayerhq.hu/ffmpeg/trunk/ ffmpeg -r 12510
echo "Done, now just do a regular configure and make to build."
echo "patching ffmpeg"
cd ffmpeg
patch -p0 <../ffmpeg.patch
cd ..
echo "linkling the dca files to ffmpeg/libavcodec"
ln -s ../../dcaenc.c ffmpeg/libavcodec/dcaenc.c
ln -s ../../dcaenc.h ffmpeg/libavcodec/dcaenc.h
echo "all done recompile ffmpeg"
