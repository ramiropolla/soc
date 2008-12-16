#! /bin/sh

echo "checking out pristine ffmpeg"
svn checkout svn://svn.mplayerhq.hu/ffmpeg/trunk/ ffmpeg -r16172

echo "patching ffmpeg"
for diff in $(ls $(pwd)/diffs/*.diff); do patch -d ffmpeg -p0 -i $diff; done

echo "copying files to libavfilter"
find -maxdepth 1 -type f -not -name $(basename $0) -exec cp {} ffmpeg/libavfilter \;
