FILES="ac3.c ac3.h ac3tab.c ac3tab.h eac3dec.c eac3.h eac3_parser.c"

echo "checking out ffmpeg svn"
for i in $FILES ac3_parser.c Makefile aac_ac3_parser.c aac_ac3_parser.h aac_parser.c allcodecs.c avcodec.h ../libavformat/allformats.h ../libavformat/raw.c; do
    rm -f ffmpeg/libavcodec/$i
done
svn checkout svn://svn.mplayerhq.hu/ffmpeg/trunk/ ffmpeg -r 9211
echo "patching ffmpeg"
cd ffmpeg
patch -p0 <../ffmpeg.patch
echo "copying the E-AC3 files to ffmpeg/libavcodec"
for i in $FILES; do
    rm -f libavcodec/$i
    ln -s ../../$i libavcodec/$i
done
echo "Done, now just do a regular configure and make to build."
