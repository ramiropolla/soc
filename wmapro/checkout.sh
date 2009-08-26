#!/bin/sh
svn co svn://svn.ffmpeg.org/ffmpeg/trunk ffmpeg
cd ffmpeg/libavcodec
rm -rf wmaprodec.c
rm -rf wmaprodata.h
ln -s ../../wmaprodec.c wmaprodec.c
ln -s ../../wmaprodata.h wmaprodata.h
cd ../
patch -p0 <../wmapro_ffmpeg.patch
patch -p0 <../audioframesize.patch

