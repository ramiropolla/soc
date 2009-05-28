#!/bin/sh
echo "checking clean ffmpeg svn"
svn checkout svn://svn.ffmpeg.org/ffmpeg/trunk ffmpeg

echo "deleting source files in libavfilter..."
rm ffmpeg/libavfilter/*.c ffmpeg/libavfilter/*.h



echo "linking SoC files where we just deleted the avfilter files from pristine"

for i in `ls *.c *.h`; do
echo "Linking $i"
ln $i ffmpeg/libavfilter/$i
done

echo "linking the SoC avfilter files into SVN..."
echo "Done, now just do a regular configure and make to build."
