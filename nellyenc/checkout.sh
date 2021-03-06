#!/bin/sh
FILES="nellymoserenc.c lowpass2.c lowpass2.h"

echo "checking out ffmpeg svn"
for i in $FILES Makefile allcodecs.c; do
    rm -f ffmpeg/libavcodec/$i
done
svn checkout svn://svn.ffmpeg.org/ffmpeg/trunk/ ffmpeg -r 14965
echo "patching ffmpeg"
cd ffmpeg
patch -p0 <../ffmpeg.patch
echo "copying the Nellymoser files to ffmpeg/libavcodec"
for i in $FILES; do
    rm -f libavcodec/$i
    ln -s ../../$i libavcodec/$i
done
echo "Done, now just do a regular configure and make to build."
