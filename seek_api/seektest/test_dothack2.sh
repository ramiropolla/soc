#!/bin/sh
echo "Test dothack.mpg from http://samples.mplayerhq.hu/MPEG2/"

./seek_test2 dothack2.mpg > res.txt

#echo "Done, please check the returned value according to the outputs."
if diff -u -w dothack2.ref res.txt ; then
    echo
    echo The output is correct.
    exit 0
else
    echo
    echo Please check the differences output.
    exit 1
fi
