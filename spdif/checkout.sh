echo "checking out ffmpeg svn"
for i in libavformat/allformats.c libavformat/Makefile libavformat/spdif.c; do
    rm -f $i
done
svn checkout svn://svn.ffmpeg.org/ffmpeg/trunk/ ffmpeg -r 19244
echo "patching ffmpeg"
cd ffmpeg
patch -p0 <../ffmpeg.patch
for i in spdif.c; do
    rm -f libavformat/$i
    ln -s ../../$i libavformat/$i
done
echo "Done, now just do a regular configure and make to build."
