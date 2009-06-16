#!/bin/sh

echo "Checking out FFmpeg SVN trunk code..."
svn co svn://svn.ffmpeg.org/ffmpeg/trunk ffmpeg &&
cd ffmpeg

echo "Adding SHA-256 support"
cp libavutil/sha1.c libavutil/sha.c
cp libavutil/sha1.h libavutil/sha.h
patch -p0 < ../sha2.patch

echo "Patching build system"
patch -p1 < ../rtmp.patch &&
cd libavformat

echo "Copying source code to libavformat dir"
cp ../../rtmp*.[ch] .

echo "Done! Now enter the ffmpeg dir, configure and make FFmpeg with some RTMP support"
