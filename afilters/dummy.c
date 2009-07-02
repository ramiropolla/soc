/* temporary development file for avfilters
* (c) 2009 Kevin DuBois <kdub432@gmail.com>
* GPL v2
*/

#include <stdio.h>
#include "avfilter.h"

typedef struct
{
    int history[100]; /*just an example */

} af_null_priv_t;


static int filter(AVFilterLink *link, AVFilterBufferRef *sample_ref)
{
    av_log(0,0, "Filter buffer\n");
    int num_samples = sample_ref->buffer->n_samples;
    int i;

    int16_t *data;
    data = (int16_t*) sample_ref->buffer->data;
    for (i=0; i < num_samples; i++)
    {
        data[i]  = data[i] +1;
    }

    return 0;
}

AVFilter avfilter_af_null =
{
    .name      = "audio_null",

    .priv_size = sizeof(af_null_priv_t),

    .inputs    = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = CODEC_TYPE_AUDIO,
                                    .filter_buffer    = filter },
                                  { .name = NULL}},

    .outputs   = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = CODEC_TYPE_AUDIO, },
                                  { .name = NULL}},
};


typedef struct
{
    int first_used;
    int last_used;

    AVFilterBufferRef buf_ref;

} af_src_priv_t;

static int av_asrc_buffer_add_samples(AVFilterContext *ctx,
        AVFilterBufferRef * samples)
{
    printf("LOAD BYTES into SRC\n");
    af_src_priv_t *priv;
    priv = ctx->priv;

    printf("last used: %i\n", priv->last_used);
    /* find the last point in the buffer */
    int first = priv->first_used;
    int last = priv->last_used;

    int attempted_load = samples->buffer->n_samples;
    if (first <= last ) /* buffer has room */
    {
        if (attempted_load > priv->buf_ref.buffer->n_samples)
        {
            /* error */
            printf("Error! not enough room\n");
            return attempted_load - priv->buf_ref.buffer->n_samples;
        }
        memcpy(&priv->buf_ref.buffer->data[last] , samples->buffer->data, attempted_load);
        priv->last_used  = priv->last_used + attempted_load;
    }


    printf("<<<<<<<<Buffer State>>>>>>>>\n");
    printf("First Used:\t%i\n",priv->first_used);
    printf("Last Used:\t%i\n", priv->last_used);

}

static int dump_next(AVFilterLink *lnk, AVFilterBufferRef *sample_ref)
{
    printf("TAKE BYTES from SRC\n");

    af_src_priv_t *priv;
    priv = (af_src_priv_t*) lnk->src->priv;


    printf("Link size is %i\n", lnk->link_size);

//    memcpy(sample_ref,  priv->buf_ref.buffer   , lnk->link_size);

    /* move samples cursor to next fresh point */
    priv->first_used = priv->first_used + lnk->link_size;

    printf("<<<<<<<<Buffer State>>>>>>>>\n");
    printf("First Used:\t%i\n",priv->first_used);
    printf("Last Used:\t%i\n", priv->last_used);

    printf("dumping buffer\n");
    return 0;

}


static int src_buf_init (AVFilterContext *ctx,
                                const char *args, void *opaque)
{
    /* allocate a fixed size for the input buffer */
    /* arbitrary value, will be modifiable */

    printf("SRC BUF INIT\n");
    af_src_priv_t *priv;
    priv = ctx->priv;

    priv->buf_ref.buffer = (AVFilterBuffer*) malloc(sizeof(AVFilterBuffer));

    priv->buf_ref.buffer->n_samples = 1024;
    priv->buf_ref.buffer->data = (int16_t *)
        calloc(priv->buf_ref.buffer->n_samples, sizeof(int16_t));

    priv->first_used = 0;
    priv->last_used =  0;

    return 0;
}



AVFilter avfilter_af_src =
{
    .name       = "audio_src",

    .priv_size  = 0,

    .init       = src_buf_init,
    .inputs     = (AVFilterPad[]) {{.name = NULL}},
    .outputs    = (AVFilterPad[]) { {.name = "default",
                                     .type = CODEC_TYPE_AUDIO,
                                     .filter_buffer = dump_next},
                                    {.name = NULL}}

};


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
    AVFilterBuffer samples;
    samples.n_samples = n_samples;
    samples.data = tbuf;
    samples.data_size = sizeof(int16_t);
    AVFilterBufferRef sample_buf;
    sample_buf.buffer = &samples;
    sample_buf.sample_type = 10;
    sample_buf.sample_rate = 128000;

    /* set up source filter */
    AVFilterContext *src_context=NULL;
    AVFilter *src_filter;
    src_filter = &avfilter_af_src;
    src_context = avfilter_open(src_filter, "filter_src");
    avfilter_register(src_filter);

    /* set up actual filter */
    AVFilterContext * avfiltcont=NULL;
    AVFilter *avfilt;
    avfilt = &avfilter_af_null;
    avfiltcont = avfilter_open(avfilt, "filterID1234");
    avfilter_register(avfilt);

    /*init filters */
    avfilter_init_filter(avfiltcont, NULL, NULL);
    avfilter_init_filter(src_context, NULL, NULL);

    /* link filters */
    avfilter_link(src_context, 0, avfiltcont, 0);
    avfilter_config_links(src_context);
    avfilter_config_links(avfiltcont);


    /* load some samples in the source filter */
    av_asrc_buffer_add_samples(src_context, &sample_buf);

    /* run this link */
    avfilter_filter_buffer(src_context->outputs[0], &sample_buf);



}

