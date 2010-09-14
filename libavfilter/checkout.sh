#! /bin/sh

echo "checking out pristine ffmpeg"
svn checkout svn://svn.ffmpeg.org/ffmpeg/trunk/ ffmpeg -r25118

echo "downloading the corresponding version of swscale"
cd ffmpeg/libswscale
svn up -r32241
cd ../..

echo "patching ffmpeg"
for diff in $(ls $(pwd)/diffs/*.diff); do patch -d ffmpeg -p0 -i $diff; done

echo "copying files to libavfilter"
cp                      \
    vf_drawbox.c        \
    vf_drawtext.c       \
    vf_fade.c           \
    vf_fps.c            \
    vf_negate.c         \
    vf_overlay.c        \
    vf_rotate.c         \
    vf_setpts.c         \
    vf_split.c          \
    vf_transpose.c      \
    vsrc_movie.c        \
    ffmpeg/libavfilter

