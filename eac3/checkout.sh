FILES="eac3dec.c ac3dec.c"

echo "checking out ffmpeg svn"
for i in $FILES ac3dec.h ac3dec_data.c ac3dec_data.h; do
    rm -f ffmpeg/libavcodec/$i
done
svn checkout svn://svn.mplayerhq.hu/ffmpeg/trunk/ ffmpeg -r 15141
echo "patching ffmpeg"
cd ffmpeg
patch -p0 <../ffmpeg.patch
echo "copying the E-AC-3 files to ffmpeg/libavcodec"
for i in $FILES; do
    rm -f libavcodec/$i
    ln -s ../../$i libavcodec/$i
done
echo "Done, now just do a regular configure and make to build."
