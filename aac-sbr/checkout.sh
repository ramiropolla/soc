#!/bin/sh

echo "Checking out FFmpeg..."
svn co svn://svn.ffmpeg.org/ffmpeg/trunk ffmpeg -r 20564
echo "Copying SBR files into FFmpeg tree..."
cp aacsbr*.[ch] ffmpeg/libavcodec
echo "Applying SBR and build system patch"
cd ffmpeg; patch -p0 < ../ffmpeg.diff
echo "All done! Now configure and make FFmpeg"
