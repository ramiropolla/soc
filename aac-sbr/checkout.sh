#!/bin/sh

echo "Checking out FFmpeg..."
svn co svn://svn.ffmpeg.org/ffmpeg/trunk ffmpeg -r 20025
echo "Copying SBR files into FFmpeg tree..."
cp aacsbr*.{c,h} ffmpeg/libavcodec
echo "Applying SBR and build system patch"
cd ffmpeg; patch -p0 < ../ffmpeg.diff
echo "All done! Now configure and make FFmpeg"
