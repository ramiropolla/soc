#!/bin/sh

echo "Checking out FFmpeg SVN trunk code..."
svn co svn://svn.mplayerhq.hu/ffmpeg/trunk ffmpeg && \
cd ffmpeg

echo "Patching build system"
patch -p1 < ../aac_enc.patch && \
cd libavcodec

echo "Copying and downloading source code to libavcodec dir"
cp ../../aacenc.c .
cp ../../aacpsy.c .
cp ../../aacpsy.h .

cd ../..
svn co svn://svn.mplayerhq.hu/soc/aac
cp aac/aactab.h ffmpeg/libavcodec/.
cp aac/aactab.c ffmpeg/libavcodec/

echo "Done! Now enter the ffmpeg dir, configure and make FFmpeg and enjoy the AAC encoder! :)"
