#!/bin/sh
echo "checking out ffmpeg svn"
svn checkout svn://svn.mplayerhq.hu/ffmpeg/trunk ffmpeg
echo "patching ffmpeg"
cd ffmpeg
patch -p0 < ../ffmpeg.patch
cd libavcodec
echo "copying the jpeg2000 files to ffmpeg/libavcodec"
ln -s ../../j2kenc.c j2kenc.c
ln -s ../../j2kdec.c j2kdec.c
ln -s ../../j2k.h j2k.h
ln -s ../../j2k.c j2k.c
ln -s ../../mqcenc.c mqcenc.c
ln -s ../../mqc.h mqc.h
ln -s ../../mqc.c mqc.c
ln -s ../../mqcdec.c mqcdec.c
ln -s ../../dwt.h dwt.h
ln -s ../../dwt.c dwt.c
echo "Done, now just do a regular configure and make to build."
