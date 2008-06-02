FILES="eac3dec.c ac3dec.c ac3dec.h ac3dec_data.c ac3dec_data.h"

echo "checking out ffmpeg svn"
for i in $FILES Makefile ac3_parser.c ac3_parser.h ac3enc.c ac3.c ac3.h; do
    rm -f ffmpeg/libavcodec/$i
done
svn checkout svn://svn.mplayerhq.hu/ffmpeg/trunk/ ffmpeg -r 13611
echo "patching ffmpeg"
cd ffmpeg
patch -p0 <../ffmpeg.patch
echo "copying the E-AC3 files to ffmpeg/libavcodec"
for i in $FILES; do
    rm -f libavcodec/$i
    ln -s ../../$i libavcodec/$i
done
echo "Done, now just do a regular configure and make to build."
