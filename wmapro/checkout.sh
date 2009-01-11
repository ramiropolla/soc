#!/bin/sh
svn co svn://svn.ffmpeg.org/ffmpeg/trunk ffmpeg
cd ffmpeg/libavcodec
ln -s ../../wma3dec.c wma3dec.c
ln -s ../../wma3data.h wma3data.h
ln -s ../../wma3.h wma3.h
cd ../
patch -p0 <../wmapro_ffmpeg.patch
patch -p0 <../audioframesize.patch

