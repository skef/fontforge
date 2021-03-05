/* Copyright (C) 2006-2012 by George Williams, 2021 by Skef Iterum */
/*
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.

 * The name of the author may not be used to endorse or promote products
 * derived from this software without specific prior written permission.

 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
 * EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 * OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <fontforge-config.h>

#include "gdraw.h"
#include "ggadgetP.h"

static GBox flowbox_box = GBOX_EMPTY; /* Don't initialize here */
static int gflowbox_inited = false;

GResInfo gflowbox_ri = {
    NULL, &ggadget_ri, NULL, NULL,
    &flowbox_box,
    NULL,
    NULL,
    NULL,
    N_("Flow Box"),
    N_("A box drawn around other gadgets to flow its contents"),
    "GGroup",
    "Gdraw",
    false,
    omf_border_type|omf_border_shape|omf_padding
    NULL,
    GBOX_EMPTY,
    NULL,
    NULL,
    NULL
};

static void _GFlowBox_Init(void) {
    if ( gflowbox_inited )
return;
    _GGadgetCopyDefaultBox(&flowbox_box);
    flowbox_box.border_type = bt_none;
    flowbox_box.border_width = 0;
    flowbox_box.padding = 0;
    _GGadgetInitDefaultBox("GFlowBox.",&flowbox_box,NULL);
    gflowbox_inited = true;
}

static void GFlowBox_destroy(GGadget *g) {
    GFlowBox *fb = (GFlowBox *) g;
    int i;

    for ( i=0; i<fb->count; ++i )
	if ( gb->children[i]!=GCD_HPad10 )
	    GGadgetDestroy(gb->children[i]);
    free(gb->children);
    _ggadget_destroy(g);
}

static void GFlowBoxMove(GGadget *g, int32 x, int32 y) {
    GFlowBox *fb = (GFlowBox *) g;
    int offx = x-g->r.x, offy = y-g->r.y;
    int i;

    for ( i=0; i<fb->count; ++i )
	if ( fb->children[i]!=GCD_HPad10 )
	    GGadgetMove(fb->children[i],
		    fb->children[i]->r.x+offx, fb->children[i]->r.y+offy);
    _ggadget_move(g,x,y);
}

struct childsizedata {
    int extra_space;		/* a default button has "extra_space" */
    int min, sized, offset;
    int opposite_min, opposite_offset;
    int rowcol;
};

struct fsizeinfo {
    struct childsizedata *childsize;
    int des, desoppo;
    int min, minoppo;
    int size, sizeoppo;
    int rowcols,
    int rowcolbaseline *;
};

static void GFlowBoxGatherMinInfo(GFlowBox *fb,struct fsizeinfo *si) {
    int c;
    GRect outer;
    int ten = GDrawPointsToPixels(fb->g.base,10);

    memset(si,0,sizeof(*si));
    si->childsize = calloc(fb->count,sizeof(struct childsizedata));

    // desired/minimum pass
    for ( c=0; c<fb->count; ++c ) {
	GGadget *g = fb->children[c];
	struct childsizedata *cs = si->childsize + c;
	if (g->state == gs_invisible)
	    continue;
	if ( g==GCD_HPad10 ) {
	    cs->min = ten;
	    cs->opposite_min = 0;
	} else {
	    GGadgetGetDesiredSize(g,&outer,NULL);
	    cs->extra_space = GBoxExtraSpace(g);
	    if (fb->vertical) {
		cs->min = outer.height;
		cs->opposite_min = outer.width;
	    } else {
		cs->min = outer.width;
		cs->opposite_min = outer.height;
	    }
	}
	si->des += cs->min;
	si->minoppo += cs->opposite_min;
	if (c+1!=fb->count) {
	    si->des += fb->vertical ? fb->vpad : fb->hpad;
	    si->minoppo += fb->vertical ? fb->hpad : fb->vpad;
	}
	if (si->desoppo < cs->opposite_min)
	    si->desoppo = cs->opposite_min;
	if (si->min < cs->min)
	    si->min = cs->min;
    }
}

static void GFlowBoxSizeTo(GFlowBox *fb,struct fsizeinfo *si, int size) {
}

static void GFlowBoxResize(GGadget *g, int32 width, int32 height) {
    struct fsizeinfo si;
    GFlowBox *fb = (GFlowBox *) g;
    int bp = GBoxBorderWidth(g->base,g->box);
    int c, size, opposize;
    int x,y;
    int old_enabled = GDrawEnableExposeRequests(g->base,false);

    GFlowBoxGatherMinInfo(fb,&si);
    width -= 2*bp; height -= 2*bp;

    size = fb->vertical ? height : width;
    if (size < si.min)
	size = si.min;

    fb->g.inner.x = fb->g.r.x + bp;
    fb->g.inner.y = fb->g.r.y + bp;

    GFlowBoxSizeTo(fb, &si, size);

    if ( si.width!=width ) {
        int vcols=0;
        for ( i=0; i<gb->cols-1; ++i )
            if ( si.cols[i].sized>0 )
                ++vcols;
        int vcols1 = vcols;
        if(si.cols[gb->cols-1].sized > 0)
            ++ vcols1;
        if ( width<si.width ) {
            for ( i=0; i<gb->cols; ++i )
                si.cols[i].sized = si.cols[i].min;
            si.width = si.minwidth;
            if ( width<si.width && gb->hpad>1 && vcols>0 ) {
                int reduce_pad = (si.width-width)/vcols + 1;
                if ( reduce_pad>gb->hpad-1 ) reduce_pad = gb->hpad-1;
                for ( i=0; i<gb->cols-1; ++i )
                    if ( si.cols[i].sized > 0 )
                        si.cols[i].sized -= reduce_pad;
                si.width -= vcols*reduce_pad;
            }
        }
        if((width > si.width) && (gb->grow_col==gb_expandglue || gb->grow_col==gb_expandgluesame )) {
            for ( i=glue_cnt=0; i<gb->cols; ++i )
                if ( si.cols[i].allglue )
                    ++glue_cnt;
            if ( glue_cnt!=0 ) {
                plus = (width-si.width)/glue_cnt;
                extra = (width-si.width-glue_cnt*plus);
                for ( i=0; i<gb->cols; ++i ) if ( si.cols[i].allglue ) {
                    si.cols[i].sized += plus + (extra>0);
                    si.width += plus + (extra>0);
                    --extra;
                }
            }
        } 
        if ((width != si.width) && gb->grow_col>=0 ) {
            int * ss = &(si.cols[gb->grow_col].sized);
            int w = si.width - *ss;
            *ss += (width-si.width);
            if(*ss < gb->hpad + 3)
                *ss = gb->hpad + 3;
            si.width = w + *ss;
        } 
        if ((width > si.width) && (vcols1!=0)) {
            plus = (width-si.width)/vcols1;
            extra = (width-si.width-vcols1*plus);
            for ( i=0; i<gb->cols; ++i ) {
                if ( si.cols[i].sized>0 ) {
                    si.cols[i].sized += plus + (extra>0);
                    si.width += plus + (extra>0);
                    --extra;
                }
            }
        }
        width = si.width;
    }

    if ( si.height!=height ) {
        int vrows=0;
        for ( i=0; i<gb->rows-1; ++i )
            if ( si.rows[i].sized>0 )
                ++vrows;
        int vrows1 = vrows;
        if(si.rows[gb->rows-1].sized > 0)
            ++ vrows1;
        if ( height<si.height ) {
            for ( i=0; i<gb->rows; ++i )
                si.rows[i].sized = si.rows[i].min;
            si.height = si.minheight;
            if ( height<si.height && gb->vpad>1 && vrows>0 ) {
                int reduce_pad = (si.height-height)/vrows + 1;
                if ( reduce_pad>gb->vpad-1 ) reduce_pad = gb->vpad-1;
                for ( i=0; i<gb->rows-1; ++i )
                    if ( si.rows[i].sized > 0 )
                        si.rows[i].sized -= reduce_pad;
                si.height -= vrows*reduce_pad;
            }
        }
        if((height > si.height) && (gb->grow_row==gb_expandglue || gb->grow_row==gb_expandgluesame )) {
            for ( i=glue_cnt=0; i<gb->rows; ++i )
                if ( si.rows[i].allglue )
                    ++glue_cnt;
            if ( glue_cnt!=0 ){
                plus = (height-si.height)/glue_cnt;
                extra = (height-si.height-glue_cnt*plus);
                for ( i=0; i<gb->rows; ++i ) if ( si.rows[i].allglue ) {
                    si.rows[i].sized += plus + (extra>0);
                    si.height += plus + (extra>0);
                    --extra;
                }
            }
        } 
        if ((height != si.height) && gb->grow_row>=0 ) {
            int * ss = &(si.rows[gb->grow_row].sized);
            int h = si.height - *ss;
            *ss += (height-si.height);
            if(*ss < gb->vpad + 3)
                *ss = gb->vpad + 3;
            si.height = h + *ss;
        } 
        if ((height > si.height) && (vrows1!=0)) {
            plus = (height-si.height)/vrows1;
            extra = (height-si.height-vrows1*plus);
            for ( i=0; i<gb->rows; ++i ) {
                if ( si.rows[i].sized>0 ) {
                    si.rows[i].sized += plus + (extra>0);
                    si.height += plus + (extra>0);
                    --extra;
                }
            }
        }
        height = si.height;
    }

    y = gb->g.inner.y;
    for ( r=0; r<gb->rows; ++r ) {
        x = gb->g.inner.x;
        for ( c=0; c<gb->cols; ++c ) {
            GGadget *g = gb->children[r*gb->cols+c];
            if ( g==GG_Glue || g==GG_ColSpan || g==GG_RowSpan || g==GCD_HPad10 )
                /* Skip it */;
            else {
                int xes, yes, es;
                totr = si.rows[r].sized;
                for ( spanr=1; r+spanr<gb->rows &&
                        gb->children[(r+spanr)*gb->cols+c]==GG_RowSpan; ++spanr )
                    totr += si.rows[r+spanr].sized;
                totc = si.cols[c].sized;
                for ( spanc=1; c+spanc<gb->cols &&
                        gb->children[r*gb->cols+c+spanc]==GG_ColSpan; ++spanc )
                    totc += si.cols[c+spanc].sized;
                if ( r+spanr!=gb->rows ) totr -= gb->vpad;
                if ( c+spanc!=gb->cols ) totc -= gb->hpad;
                es = GBoxExtraSpace(g);
                xes = si.cols[c].extra_space - es;
                yes = si.rows[r].extra_space - es;
                if ( g->state!=gs_invisible )
                    GGadgetResize(g,totc-2*xes,totr-2*yes);
                GGadgetMove(g,x+xes,y+yes);
            }
            x += si.cols[c].sized;
        }
        y += si.rows[r].sized;
    }

    free(si.cols); free(si.rows);

    gb->g.inner.width = width; gb->g.inner.height = height;
    gb->g.r.width = width + 2*bp; gb->g.r.height = height + 2*bp;
    GDrawEnableExposeRequests(g->base,old_enabled);
    GDrawRequestExpose(g->base,&g->r,false);
}

static void GFlowBoxGetDesiredSize(GGadget *g, GRect *outer, GRect *inner) {
    struct fsizeinfo si;
    GFlowBox *fb = (GFlowBox *) g;
    int bp = GBoxBorderWidth(g->base,g->box);
    int width, height;

    GFlowBoxGatherMinInfo(fb,&si);
    if ( fb->vertical ) {
	if ( g->desired_height>0 ) {
	    GFlowBoxSizeTo(fb, &si, g->desired_height);
	    height = si.size;
	    width = si.sizeoppo;
	} else {
	    height = si.des;
	    width = si.desoppo;
	}
    } else {
	if ( g->desired_width>0 ) {
	    GFlowBoxSizeTo(fb, &si, g->desired_width);
	    width = si.size;
	    height = si.sizeoppo;
	} else {
	    width = si.des;
	    height = si.desoppo;
	}
    }

    if ( inner!=NULL ) {
	inner->x = inner->y = 0;
	inner->width = width; inner->height = height;
    }
    if ( outer!=NULL ) {
	outer->x = outer->y = 0;
	outer->width = width+2*bp; outer->height = height+2*bp;
    }
    free(si.childsize); free(si.rowcolbaseline);
}

static int GFlowBoxFillsWindow(GGadget *g) {
return( true );
}

static int expose_nothing(GWindow pixmap, GGadget *g, GEvent *event) {
    if ( g->state!=gs_invisible )
	GBoxDrawBorder(pixmap,&g->r,g->box,g->state,false);
    return true;
}

struct gfuncs gflowbox_funcs = {
    0,
    sizeof(struct gfuncs),

    expose_nothing,	/* Expose */
    _ggadget_noop,	/* Mouse */
    _ggadget_noop,	/* Key */
    NULL,
    NULL,		/* Focus */
    NULL,
    NULL,

    _ggadget_redraw,
    GFlowBoxMove,
    GFlowBoxResize,
    _ggadget_setvisible,
    _ggadget_setenabled,
    _ggadget_getsize,
    _ggadget_getinnersize,

    GFlowBox_destroy,

    NULL,
    NULL,
    NULL,
    NULL,
    NULL,

    NULL,
    NULL,

    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,

    GFlowBoxGetDesiredSize,
    _ggadget_setDesiredSize,
    GFlowBoxFillsWindow,
    NULL
};

/*void GHVBoxSetExpandableCol(GGadget *g,int col) {
    GHVBox *gb = (GHVBox *) g;
    if ( col<gb->cols )
	gb->grow_col = col;
}*/

GGadget *GFlowBoxCreate(struct gwindow *base, GGadgetData *gd, void *data) {
    GFlowBox *fb = calloc(1,sizeof(GFlowBox));

    if ( !gflowbox_inited )
	_GFlowBox_Init();

    for ( fb->count=0; gd->u.boxelements[fb->count]!=NULL; ++fb->count );

    fb->g.funcs = &gflowbox_funcs;
    _GGadget_Create(&fb->g,base,gd,data,&flowbox_box);
    fb->hpad = fb->vpad = GDrawPointsToPixels(base,2);

    fb->g.takes_input = false;
    fb->g.takes_keyboard = false;
    fb->g.focusable = false;

    fb->children = malloc(fb->count*sizeof(GGadget *));
    for ( int i=0; i<fb->count; ++i) {
	GGadgetCreateData *gcd = gd->u.boxelements[i];
	if ( gcd==GCD_HPad10 )
	    gb->children[i] = (GGadget *) gcd;
	else {
	    gcd->gd.pos.x = gcd->gd.pos.y = 1;
	    fb->children[i] = gcd->ret = (gcd->creator)(base,&gcd->gd,gcd->data);
	    gcd->ret->contained = true;
	}
    }
    return fb;
}

GResInfo *_GFlowBoxRIHead(void) {
    _GFlowBox_Init();
    return &gflowbox_ri;
}
