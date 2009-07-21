#include <stdio.h>
#include "avfilter.h"

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
        memcpy(&priv->buf_ref.buffer->data[last] , samples->buffer->data,
                sizeof(int16_t) * attempted_load);
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


static int query_af_src_formats(AVFilterContext *ctx)
{
    av_log(0,0, "query formats\n");
    AVFilterFormats *formats1;
    formats1 = avfilter_all_sampleformats();
    avfilter_set_common_formats(ctx,formats1);
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
                                    {.name = NULL}},

    .query_formats = query_af_src_formats
};

