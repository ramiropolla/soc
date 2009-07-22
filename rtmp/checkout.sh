#!/bin/sh

echo "Checking out FFmpeg SVN trunk code..."
svn co svn://svn.ffmpeg.org/ffmpeg/trunk ffmpeg &&
cd ffmpeg

echo "Patching build system"
patch -p1 < ../rtmp.patch &&
cd libavformat

echo "Copying source code to libavformat dir"
cp ../../rtmp*.[ch] .

echo "Done! Now enter the ffmpeg dir, configure and make FFmpeg with some RTMP support"
