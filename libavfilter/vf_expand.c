#include "avfilter.h"

typedef struct
{
    unsigned expw, exph;
} priv_t;

AVFrame *get_buffer(AVFilterContext *ctx, unsigned w, unsigned h, int fmt)
{
    priv_t *p = ctx->priv;
    AVFrame *buf;

    if((buf = ctx->out->dst->filt->get_buffer(p->expw, p->exph, fmt))) {
    }

    return buf;
}

AVFilter vf_expand =
{
    .name = "expand",
    .author = "Bobby Bingham",
    .get_buffer = get_buffer,
};
