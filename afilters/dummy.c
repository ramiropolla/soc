/* temporary development file for avfilters
* (c) 2009 Kevin DuBois <kdub432@gmail.com>
*/

#include "af_null.c" /*FIXME: bad, i know. prototyping :) */


int main()
{

    printf("AVfilter version: %x\n", avfilter_version());


    /* Simulates a 1024 buffer of sl16 audio data */
    /* temporary setup, probably a cleaner way i want to do all this */
    int16_t * tbuf;
    int i, n_samples = 1024;
    tbuf = calloc(n_samples, sizeof(int16_t));
    for(i=0;i<n_samples;i++)
    {
        tbuf[i] = (int16_t) 100 * sin(2*3.141/100 * i);
    }   // sine wave, period 1024/100, range, 100 -> -100

    AVFilterSamples samples;
    samples.n_samples = n_samples;
    samples.data = tbuf;
    samples.data_size = sizeof(int16_t);

    AVFilterSamplesRef sample_buf;
    sample_buf.samples = &samples;
    sample_buf.sample_type = 10;
    sample_buf.sample_rate = 128000;

    /* avfilter context */
    AVFilterContext * avfiltcont;
    AVFilter *avfilt;

    /*set up avfilter */
    avfilt = &avfilter_af_null;

    /* this should initialize the avfiltcont */
    printf("Opening avfilter\n");
    avfiltcont = avfilter_open(avfilt, "kevinfilter");
    printf("avfilter_open done\n");


    /* Register filters*/
    avfilter_register(avfilt);

    /* initialize the filters */
    printf("Starting filter chain\n");
    avfilter_init_filter(avfiltcont, NULL, NULL);

    printf("Alright, we got %i inputs and %i outputs\n",
            avfiltcont->input_count,
            avfiltcont->output_count);

    /* run the filter */
    printf("Running filter chain\n");
    /* FIXME: trying to run segfaults :(. didnt set up linking right, i think*/
    avfilter_start_buffer(avfiltcont->inputs, &sample_buf);
    //avfilter_start_buffer(AVFilterLink *link, AVFilterPicRef *picref);



    /* uninit things */

}

