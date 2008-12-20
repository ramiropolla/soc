echo "checking out ffmpeg svn"
svn checkout svn://svn.ffmpeg.org/ffmpeg/trunk/ ffmpeg
echo "copying the ac3 decoder files to ffmpeg/libavcodec"
cp ac3_decoder.c ffmpeg/libavcodec/
cp ac3_decoder.h ffmpeg/libavcodec/
echo "patching ffmpeg"
cd ffmpeg ; patch -p0 <../lgpl_ac3_decoder.patch ; cd ..
echo "Done, now just do a regular configure and make to build."
