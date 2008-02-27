/*
 * Movie file video source filter
 * Copyright (c) 2008 Victor Paesa
 *
 * This file is part of FFmpeg.
 *
 * FFmpeg is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * FFmpeg is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with FFmpeg; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

 /*
 # One usage example follows, usable too as test scenario.
 #TODO Eventually move it into FFmpeg docs.

;----------------------------movie.desc------------------------------
;
; Parameters of movie filter are
; seekpoint in microseconds : string format : string filename
;
; We can overlay a second movie on top of a main one
;
; input -----------> deltapts0 --> overlay --> output
;                                    ^
; movie --> scale--> deltapts1 ------|
;

[filters]
movie=movie=3200000:avi:input.avi

; We make both the main movie and the overlaid one start at the same time
deltapts0=setpts=PTS-STARTPTS
deltapts1=setpts=PTS-STARTPTS

scale=scale=180:144
overlay=overlay=16:16

[links]
deltapts0:0=overlay:0

movie:0=scale:0
scale:0=deltapts1:0
deltapts1:0=overlay:1

[inputs]
default=deltapts0:0

[outputs]
default=overlay:0
;----------------------------movie.desc------------------------------

# Then launch it from command line
ffmpeg -i input.avi -vfilters graph_file=movie.desc output.avi

 */

#include "avformat.h"
#include "avfilter.h"

typedef struct {
    // Filter parameters
    int64_t           seek_point; //< seekpoint in microseconds
    char              format_name[16];
    char              file_name[255];
    // Needed to load movies
    AVFormatContext  *pFormatCtx;
    int               videoStream;
    AVCodecContext   *pCodecCtx;
    int               is_done;
    AVFrame          *pFrame;

    int w, h;
    AVFilterPicRef *pic;
} MovieContext;

int movie_init(AVFilterContext *ctx)
{
    AVInputFormat  *file_iformat = NULL;
    int             i;
    AVCodec        *pCodec;
    int64_t         timestamp;
    MovieContext   *mv = ctx->priv;

    /* av_log(ctx, AV_LOG_INFO, "movie_init() seek:%lld format:%s file:%s\n",
        mv->seek_point, mv->format_name, mv->file_name); */
    av_register_all();
    // Try to find the movie format (container)
    if (*mv->format_name)
        file_iformat = av_find_input_format(mv->format_name);
    else
        file_iformat = NULL;
    mv->pFormatCtx = NULL;
    if (av_open_input_file(&mv->pFormatCtx, mv->file_name, file_iformat, 0, NULL) != 0) {
        av_log(ctx, AV_LOG_ERROR,
            "movie_init() Failed to av_open_input_file '%s'\n", mv->file_name);
        return -1;
    }
    if(av_find_stream_info(mv->pFormatCtx)<0) {
        av_log(ctx, AV_LOG_ERROR, "movie_init() Failed to find stream info\n");
        return -1;
    }

    // if seeking requested, we execute it
    if (mv->seek_point > 0) {
        timestamp = mv->seek_point;
        // add the stream start time, should it exist
        if (mv->pFormatCtx->start_time != AV_NOPTS_VALUE)
            timestamp += mv->pFormatCtx->start_time;
        if (av_seek_frame(mv->pFormatCtx, -1, timestamp, AVSEEK_FLAG_BACKWARD) < 0) {
            av_log(ctx, AV_LOG_ERROR, "%s: could not seek to position %"PRId64"\n",
                mv->file_name, timestamp);
        }
    }

    // To make things nice and easy, we simply use the first video stream we find
    // TODO: allow to choose the video stream
    mv->videoStream = -1;
    for(i = 0; i < mv->pFormatCtx->nb_streams; i++)
        if(mv->pFormatCtx->streams[i]->codec->codec_type==CODEC_TYPE_VIDEO) {
            mv->videoStream = i;
            break;
        }
    if(mv->videoStream == -1) {
        av_log(ctx, AV_LOG_ERROR, "movie_init() No video stream found\n");
        return -1;
    }
    // Get a pointer to the codec context for the video stream
    mv->pCodecCtx = mv->pFormatCtx->streams[mv->videoStream]->codec;

    /*
     * So now we've got a pointer to the so-called codec context for our video
     * stream, but we still have to find the actual codec and open it.
     */
    // Find the decoder for the video stream
    pCodec = avcodec_find_decoder(mv->pCodecCtx->codec_id);
    if(!pCodec) {
        av_log(ctx, AV_LOG_ERROR, "movie_init() Failed to find any codec\n");
        return -1;
    }

    // Open codec
    if(avcodec_open(mv->pCodecCtx, pCodec)<0) {
        av_log(ctx, AV_LOG_ERROR, "movie_init() Failed to open codec\n");
        return -1;
    }

    // Allocate a video frame to store the decoded images in.
    if(! (mv->pFrame = avcodec_alloc_frame()) ) {
        av_log(ctx, AV_LOG_ERROR, "movie_init() Failed to alloc frame\n");
        return -1;
    }

    mv->w = mv->pCodecCtx->width;
    mv->h = mv->pCodecCtx->height;

    return 0;
}

static int init(AVFilterContext *ctx, const char *args, void *opaque)
{
    MovieContext *mv = ctx->priv;

    if(args) {
        int num_fields = sscanf(args, "%"PRId64":%15[^:]:%255s",
                            &mv->seek_point, mv->format_name, mv->file_name);
        if (3 == num_fields)
            /* av_log(ctx, AV_LOG_INFO,
                "init() args:'%s'\n\tseek:%lld format:%s name:%s\n",
                args, mv->seek_point, mv->format_name, mv->file_name); */
            // sanity check parms
            if (mv->seek_point >= 0 && *mv->file_name)
                return movie_init(ctx);
        }
        else
            av_log(ctx, AV_LOG_ERROR, "init() expected 3 arguments:'%s'\n", args);
    return -1;
}

static int query_formats(AVFilterContext *ctx)
{
    MovieContext *mv = ctx->priv;

    avfilter_set_common_formats(ctx,
        avfilter_make_format_list(1, mv->pCodecCtx->pix_fmt));
    return 0;
}

static int config_props(AVFilterLink *link)
{
    MovieContext *mv = link->src->priv;

    link->w = mv->w;
    link->h = mv->h;

    return 0;
}

int movie_get_frame(AVFilterLink *link)
{
    AVPacket packet;
    int      frameFinished;

    MovieContext *mv = link->src->priv;

    if (mv->is_done == 1) return 0;

    if(!mv->pic)
        mv->pic = avfilter_get_video_buffer(link, AV_PERM_WRITE );
    //av_log(link->src, AV_LOG_INFO, "movie_get_frame() w:%d h:%d\n", mv->w, mv->h);

    // Get frame
    while(av_read_frame(mv->pFormatCtx, &packet)>=0)
    {
        // Is this a packet from the video stream?
        if(packet.stream_index == mv->videoStream)
        {
            // Decode video frame
            avcodec_decode_video(mv->pCodecCtx, mv->pFrame, &frameFinished,
                packet.data, packet.size);

            // Did we get a video frame?
            if(frameFinished)
            {
                memcpy(mv->pic->data,     mv->pFrame->data,
                       sizeof(mv->pFrame->data));
                memcpy(mv->pic->linesize, mv->pFrame->linesize,
                       sizeof(mv->pFrame->linesize));

                // Advance in the time line
                mv->pic->pts = av_rescale_q(packet.pts,
                    mv->pFormatCtx->streams[mv->videoStream]->time_base,
                    AV_TIME_BASE_Q);
                /* av_log(link->src, AV_LOG_INFO,
                  "movie_get_frame(%s) packet pts:%lld %lf vfpts:%lld\n",
                  mv->file_name, packet.pts, (double)packet.pts *
                  av_q2d(mv->pFormatCtx->streams[mv->videoStream]->time_base),
                  mv->pic->pts);*/

                // We got it. Free the packet since we are returning
                av_free_packet(&packet);

                return 0;
            }
        }
        // Free the packet that was allocated by av_read_frame
        av_free_packet(&packet);
    }
    // On multi-frame source we should stop the mixing process when
    // the movie source does not have more frames
    mv->is_done = 1;
    return 0;
}

static int request_frame(AVFilterLink *link)
{
    AVFilterPicRef *out;
    MovieContext *mv = link->src->priv;

    movie_get_frame(link);

    out = avfilter_ref_pic(mv->pic, ~AV_PERM_WRITE);
    out->pixel_aspect = mv->pCodecCtx->sample_aspect_ratio;

    avfilter_start_frame(link, out);
    avfilter_draw_slice(link, 0, out->h);
    avfilter_end_frame(link);

    return 0;
}

static void uninit(AVFilterContext *ctx)
{
    MovieContext *mv = ctx->priv;

    if(mv->pic)
        avfilter_unref_pic(mv->pic);
}

AVFilter avfilter_vsrc_movie =
{
    .name      = "movie",
    .priv_size = sizeof(MovieContext),
    .query_formats = query_formats,

    .init      = init,
    .uninit    = uninit,

    .inputs    = (AVFilterPad[]) {{ .name = NULL }},
    .outputs   = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = CODEC_TYPE_VIDEO,
                                    .request_frame   = request_frame,
                                    .config_props    = config_props, },
                                  { .name = NULL}},
};

