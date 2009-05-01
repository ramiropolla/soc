#! /bin/sh

echo "checking out pristine ffmpeg"
svn checkout svn://svn.ffmpeg.org/ffmpeg/trunk/ ffmpeg -r18726

echo "patching ffmpeg"
for diff in $(ls $(pwd)/diffs/*.diff); do patch -d ffmpeg -p0 -i $diff; done

echo "copying files to libavfilter"
find $(pwd) -maxdepth 1 -type f -not -name $(basename $0) -exec cp {} ffmpeg/libavfilter \;

echo "copying libavfilter doc files"
cp -p $(pwd)/doc/* ffmpeg/doc
