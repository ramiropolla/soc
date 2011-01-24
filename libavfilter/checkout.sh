#! /bin/sh

echo "checking out pristine ffmpeg"
svn checkout svn://svn.ffmpeg.org/ffmpeg/trunk@26400 ffmpeg

echo "downloading the corresponding version of swscale"
svn checkout svn://svn.ffmpeg.org/mplayer/trunk/libswscale@32676 ffmpeg/libswscale

echo "patching ffmpeg"
for diff in $(ls $(pwd)/diffs/*.diff); do patch -d ffmpeg -p0 -i $diff; done

echo "copying files to libavfilter"
cp                      \
    vf_drawtext.c       \
    vf_fade.c           \
    vf_fps.c            \
    vf_negate.c         \
    vf_rotate.c         \
    vf_split.c          \
    vsrc_movie.c        \
    ffmpeg/libavfilter

