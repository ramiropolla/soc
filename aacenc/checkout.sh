#!/bin/sh

echo "Checking out FFmpeg SVN trunk code..."
svn co svn://svn.ffmpeg.org/ffmpeg/trunk ffmpeg &&
cd ffmpeg

echo "Patching build system"
patch -p1 < ../aac_enc.patch &&
cd libavcodec

echo "Copying and downloading source code to libavcodec dir"
cp ../../aac.h .
cp ../../aacenc.c .
cp ../../aacenc.h .
cp ../../aaccoder.c .
cp ../../aacpsy.c .
cp ../../psymodel.c .
cp ../../psymodel.h .

echo "Done! Now enter the ffmpeg dir, configure and make FFmpeg and enjoy the AAC encoder! :)"
