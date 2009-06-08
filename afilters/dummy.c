/* temporary development file for avfilters
* (c) 2009 Kevin DuBois <kdub432@gmail.com>
*/

#include "avfilter.h"
#include "af_null.c"

int main()
{

    printf("%x\n", avfilter_version());


    /* Simulates a 1024 buffer of sl16 audio data */
    int16_t * tbuf;
    tbuf = calloc(1024, sizeof(int16_t));
    int i;
    for(i=0;i<1024;i++)
    {
        tbuf[i] = (int16_t) 100 * sin(2*3.141/100 * i);
        printf("%i\n", tbuf[i]);

    }   // sine wave, period 1024/100, range, 100 -> -100


    /* avfilter context */
    AVFilterContext * avfiltcont;

    AVFilter *avfilt=NULL, *filt2=NULL;


    avfilt = &avfilter_af_null;


    AVFilterPad *pad;
    avfilter_insert_inpad(0, avfiltcont, pad );

    avfiltcont = avfilter_open(avfilt, "kevinfilter");

    avfilter_register(avfilt);
    filt2 = avfilter_get_by_name("kevinfilter");

    /*run buffer now*/



}

