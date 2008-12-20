#!/bin/sh
svn co svn://svn.ffmpeg.org/ffmpeg/trunk ffmpeg
cd ffmpeg/libavcodec
ln -s ../../wma3dec.c wma3dec.c
cd ../
patch -p0 <../wmapro_ffmpeg.patch

