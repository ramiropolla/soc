/* temporary development file for avfilters
* (c) 2009 Kevin DuBois <kdub432@gmail.com>
* GPL v2
*/

#include <stdio.h>
#include "avfilter.h"
#include "af_src.c"
#include "af_null.c"
#include "af_vol.c"

#if 0
int dump_avfiltlink(AVFilterLink *link)
{
    if (!link)
        return 0;

    printf("\tLink dump...\n");
    printf("\tSource:\t0x%x\n", link->src);
    printf("\tDest:\t0x%x\n", link->dst);


    switch (link->init_state)
    {
    case AVLINK_UNINIT:
        printf("\tState: AVLINK_UNINIT\n");
        break;
    case AVLINK_STARTINIT:
        printf("\tState: AVLINK_STARTINIT\n");
        break;
    case AVLINK_INIT:
        printf("\tState: AVLINK_INIT\n");
        break;
    default:
        printf("\tState: ???\n");
        break;
    }

}

int dump_avfiltcont(AVFilterContext *cont)
{
    printf("\n--------------AVFILTCONT DUMP-------------\n");
    if (!cont)
    {
        printf("Error, null argument\n");
        printf("------------END AVFILTCONT DUMP-----------\n\n");
        return -1;
    }

    printf("Cont addr:\t%X\n", cont);
    printf("Class:\t\t%x\n", cont->av_class);
    printf("Filter:\t\t%x\n", cont->filter);
    printf("Name:\t\t%s\n", cont->name);
    printf("Input Count:\t%i\n", cont->input_count);
    printf("Input Pads:\t%x\n", cont->input_pads);
    printf("Input Links:\t%x\n", cont->inputs);
//    dump_avfiltlink(cont->inputs[0]);
    printf("Output Count:\t%i\n", cont->output_count);
    printf("Output Pads:\t%x\n", cont->output_pads);
    printf("Output Links:\t%x\n", cont->outputs);
    dump_avfiltlink(cont->outputs[0]);

    printf("------------END AVFILTCONT DUMP-----------\n\n");
    return 0;
}
#endif

int main()
{
    printf("AVfilter version: %x\n", avfilter_version());

    /* Simulates a 1024 buffer of sl16 audio data */
    /* temporary setup, probably a cleaner way i want to do all this */
    int16_t * tbuf;
    int i, n_samples = 512;
    tbuf = calloc(n_samples, sizeof(int16_t));
    for(i=0;i<n_samples;i++)
    {
#define SINEWAVE 0
#if SINEWAVE
        tbuf[i] = (int16_t) 100 * sin(2*3.141/100 * i);
#else
        tbuf[i] = i;
#endif
    }   // sine wave, period 1024/100, range, 100 -> -100

    /* set up formats */

    AVFilterBuffer samples;
    AVFilterBufferRef sample_buf;
    samples.n_samples = n_samples;
    samples.data = tbuf;
    samples.data_size = sizeof(int16_t);

    sample_buf.sample_type = 10;
    sample_buf.buffer = &samples;
    sample_buf.sample_rate = 128000;

    /* set up source filter */
    AVFilterContext *src_context=NULL;
    AVFilter *src_filter;
    src_filter = &avfilter_af_src;
    src_context = avfilter_open(src_filter, "filter_src");
    avfilter_register(src_filter);

    /* set up the first filter */
    AVFilterContext * avfiltcont=NULL;
    AVFilter *avfilt;
    avfilt = (AVFilter*)malloc(sizeof(AVFilter));
    memcpy(avfilt, &avfilter_af_null, sizeof(AVFilter));
    avfiltcont = avfilter_open(avfilt, "filterID1234");
    avfilter_register(avfilt);

    /* set up the first filter */
    AVFilterContext * avfiltcont2=NULL;
    AVFilter *avfilt2;
    avfilt2 = (AVFilter*)malloc(sizeof(AVFilter));
    memcpy(avfilt2, &avfilter_af_volume, sizeof(AVFilter));
    avfiltcont2 = avfilter_open(avfilt2, "filtery");
    avfilter_register(avfilt2);

    /*init filters */
    avfilter_init_filter(avfiltcont, NULL, NULL);
    avfilter_init_filter(avfiltcont2, NULL, NULL);
    avfilter_init_filter(src_context, NULL, NULL);

    /* configure formats?  */



    /* merge lists? */
    AVFilterFormats * merge;

//    merge = avfilter_merge_formats(AVFilterFormats *a, AVFilterFormats *b)


    /* link filters */
    avfilter_link(src_context, 0, avfiltcont, 0);
    avfilter_link(avfiltcont, 0, avfiltcont2, 0);
    avfilter_config_links(src_context);
    avfilter_config_links(avfiltcont);
    avfilter_config_links(avfiltcont2);





    /*dump_avfiltcont(src_context);
    dump_avfiltcont(avfiltcont);
    dump_avfiltcont(avfiltcont2);*/

    /* load some samples in the source filter */
    av_asrc_buffer_add_samples(src_context, &sample_buf);

    /* run this link */
    printf("running first link\n");
    avfilter_filter_buffer(src_context->outputs[0], &sample_buf);
    printf("running second link\n");
    avfilter_filter_buffer(avfiltcont->outputs[0], &sample_buf);

}

