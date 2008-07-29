#!/bin/sh
for file in `ls *.wv`
do
    for br in 32k 48k 64k 96k 128k 160k 192k
    do
        ffmpeg -i $file -ab $br -y "`basename $file .wv`.$br.mp4"
    done
done
