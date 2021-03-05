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

#define GG_Pad10       ((GGadget *) -4)		// Must match GCD_Hpad10 in ggadget.h

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
    omf_border_type|omf_border_shape|omf_padding,
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
    int c;

    for ( c=0; c<fb->count; ++c )
	if ( fb->children[c]!=GG_Pad10 )
	    GGadgetDestroy(fb->children[c]);
    free(fb->children);
    _ggadget_destroy(g);
}

static void GFlowBoxMove(GGadget *g, int32 x, int32 y) {
    GFlowBox *fb = (GFlowBox *) g;
    int offx = x-g->r.x, offy = y-g->r.y;
    int c;

    for ( c=0; c<fb->count; ++c )
	if ( fb->children[c]!=GG_Pad10 )
	    GGadgetMove(fb->children[c],
		    fb->children[c]->r.x+offx, fb->children[c]->r.y+offy);
    _ggadget_move(g,x,y);
}

struct childsizedata {
    int extra_space;		/* a default button has "extra_space" */
    int min, size, offset;
    int oppomin;
    int opposize;
    int oppooffset; // Relative to rowcol
    int rowcol;
};

struct fsizeinfo {
    struct childsizedata *childsize;
    int *rowcolbaseline;
    int des, oppodes;
    int min, oppomin;
    int size, opposize;
    int rowcols;
};

static void GFlowBoxGatherMinInfo(GFlowBox *fb,struct fsizeinfo *si) {
    GRect outer;
    int c, ten = GDrawPointsToPixels(fb->g.base,10);

    memset(si,0,sizeof(*si));
    si->childsize = calloc(fb->count,sizeof(struct childsizedata));

    // desired/minimum pass
    for ( c=0; c<fb->count; ++c ) {
	GGadget *g = fb->children[c];
	struct childsizedata *cs = si->childsize + c;
	if (g->state == gs_invisible)
	    continue;
	if ( g==GG_Pad10 ) {
	    cs->min = ten;
	    cs->oppomin = 0;
	} else {
	    GGadgetGetDesiredSize(g,&outer,NULL);
	    cs->extra_space = GBoxExtraSpace(g);
	    if (fb->vertical) {
		cs->min = outer.height;
		cs->oppomin = outer.width;
	    } else {
		cs->min = outer.width;
		cs->oppomin = outer.height;
	    }
	}
	si->des += cs->min;
	si->oppomin += cs->oppomin;
	if (c+1!=fb->count) {
	    si->des += fb->vertical ? fb->vpad : fb->hpad;
	    si->oppomin += fb->vertical ? fb->hpad : fb->vpad;
	}
	if (si->oppodes < cs->oppomin)
	    si->oppodes = cs->oppomin;
	if (si->min < cs->min)
	    si->min = cs->min;
    }
}

static void GFlowBoxSizeTo(GFlowBox *fb,struct fsizeinfo *si, int size) {
    int rowcolstart=0, rowcoloffset=0, cursize=0, opposize=0;
    int pad = fb->vertical ? fb->vpad : fb->hpad;
    int oppopad = fb->vertical ? fb->hpad : fb->vpad;
    int c, j;

    si->rowcols = 0;
    if (si->rowcolbaseline!=NULL)
	free(si->rowcolbaseline);
    si->rowcolbaseline = calloc(fb->count,sizeof(int));

    for (c=0; c<fb->count; ++c) {
	struct childsizedata *cs = si->childsize + c;
	if (c!=rowcolstart && cursize + cs->min + pad>size) {
	    for (j=rowcolstart; j<c; ++j) {
		;
	    }
	    si->rowcolbaseline[si->rowcols] = rowcoloffset;
	    si->rowcols++;
	    rowcoloffset += opposize + oppopad;
	    opposize = 0;
	    rowcolstart = c;
	    cursize = 0;
	}
	if (opposize < cs->oppomin)
	    opposize = cs->oppomin;
	if (c!=rowcolstart)
	    cursize += pad;
	cs->offset = cursize;
	cs->size = si->min; // XXX
	cs->opposize = si->oppomin; // XXX
	cs->rowcol = si->rowcols;
	cursize += cs->min;
    }
    for (j=rowcolstart; j<c; ++j) {
	;
    }
    si->rowcolbaseline[si->rowcols] = rowcoloffset;
    si->rowcols++;
    si->size = size;
    si->opposize = rowcoloffset + opposize;
}

static void GFlowBoxResize(GGadget *g, int32 width, int32 height) {
    struct fsizeinfo si;
    GFlowBox *fb = (GFlowBox *) g;
    int bp = GBoxBorderWidth(g->base,g->box);
    int c, size;
    int x,y,cw,ch;
    int old_enabled = GDrawEnableExposeRequests(g->base,false);

    GFlowBoxGatherMinInfo(fb,&si);
    width -= 2*bp; height -= 2*bp;

    size = fb->vertical ? height : width;
    if (size < si.min)
	size = si.min;

    fb->g.inner.x = fb->g.r.x + bp;
    fb->g.inner.y = fb->g.r.y + bp;

    GFlowBoxSizeTo(fb, &si, size);

    for ( c=0; c<fb->count; ++c ) {
	GGadget *g = fb->children[c];
	struct childsizedata *cs = si.childsize + c;
	if (g==GG_Pad10)
	    continue;
	if ( fb->vertical ) {
	    x = cs->oppooffset + si.rowcolbaseline[cs->rowcol];
	    y = cs->offset;
	    cw = cs->opposize;
	    ch = cs->size;
	} else {
	    x = cs->offset;
	    y = cs->oppooffset + si.rowcolbaseline[cs->rowcol];
	    cw = cs->size;
	    ch = cs->opposize;
	}
	if ( g->state!=gs_invisible )
	    GGadgetResize(g,cw,ch);
	GGadgetMove(g,fb->g.inner.x+x,fb->g.inner.y+y);
    }

    free(si.childsize); free(si.rowcolbaseline);

    fb->g.inner.width = fb->vertical ? si.opposize : si.size;
    fb->g.inner.height = fb->vertical ? si.size : si.opposize;
    fb->g.r.width = fb->g.inner.width + 2*bp;
    fb->g.r.height = fb->g.inner.height + 2*bp;
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
	    width = si.opposize;
	} else {
	    height = si.des;
	    width = si.oppodes;
	}
    } else {
	if ( g->desired_width>0 ) {
	    GFlowBoxSizeTo(fb, &si, g->desired_width);
	    width = si.size;
	    height = si.opposize;
	} else {
	    width = si.des;
	    height = si.oppodes;
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
	    fb->children[i] = (GGadget *) gcd;
	else {
	    gcd->gd.pos.x = gcd->gd.pos.y = 1;
	    fb->children[i] = gcd->ret = (gcd->creator)(base,&gcd->gd,gcd->data);
	    gcd->ret->contained = true;
	}
    }
    return &fb->g;
}

GResInfo *_GFlowBoxRIHead(void) {
    _GFlowBox_Init();
    return &gflowbox_ri;
}
