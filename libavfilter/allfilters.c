/*
 * filter registration
 * copyright (c) 2008 Vitor Sessak
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

#include "avfilter.h"

#define REGISTER_VF(X,x) { \
          extern AVFilter avfilter_vf_##x ; \
          if(ENABLE_VF_##X )  avfilter_register(&avfilter_vf_##x ); }


#define REGISTER_VSRC(X,x) { \
          extern AVFilter avfilter_vsrc_##x ; \
          if(ENABLE_VSRC_##X )  avfilter_register(&avfilter_vsrc_##x ); }

void avfilter_register_all(void)
{
    static int initialized;

    if (initialized)
        return;
    initialized = 1;

    REGISTER_VF(CROP,crop);
    REGISTER_VF(DRAWBOX,drawbox);
    REGISTER_VF(FIFO,fifo);
    REGISTER_VF(FORMAT,format);
    REGISTER_VF(FPS,fps);
    REGISTER_VF(GRAPH,graph);
    REGISTER_VF(GRAPHDESC,graphdesc);
    REGISTER_VF(GRAPHFILE,graphfile);
    REGISTER_VF(HFLIP,hflip);
    REGISTER_VF(NEGATE,negate);
    REGISTER_VF(NOFORMAT,noformat);
    REGISTER_VF(OVERLAY,overlay);
    REGISTER_VF(ROTATE,rotate);
    REGISTER_VF(SCALE,scale);
    REGISTER_VF(SETPTS,setpts);
    REGISTER_VF(SLICIFY,slicify);
    REGISTER_VF(SPLIT,split);
    REGISTER_VF(TRANSPOSE,transpose);
    REGISTER_VF(VFLIP,vflip);

    REGISTER_VSRC(MOVIE,movie);
}
