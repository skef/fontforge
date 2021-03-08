/* Copyright (C) 2000-2012 by George Williams */
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
#include "gkeysym.h"
#include "gwidget.h"
#include "ustring.h"

static int GScroll1BoxTime = 500;		/* half a second between scrolls when mouse out of listbox */

static GBox scroll1box_box = GBOX_EMPTY; /* Don't initialize here */;
static int gscroll1box_inited = false;

static GResInfo gscroll1box_ri = {
    NULL, &ggadget_ri,NULL, NULL,
    &scroll1box_box,
    NULL,
    NULL,
    NULL,
    N_("Scroll Box"),
    N_("A box used to scroll widgets vertically or horizontally"),
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

static void GScroll1Box_Init() {
    if ( gscroll1box_inited )
	return;
    _GGadgetCopyDefaultBox(&scroll1box_box);
    scroll1box_box.border_type = bt_none;
    scroll1box_box.border_width = 0;
    scroll1box_box.padding = 0;
    _GGadgetInitDefaultBox("GScroll1Box",&scroll1box_box, NULL);
    gscroll1box_inited = true;
}

static void GScroll1BoxDestroy(GGadget *g) {
    GScroll1Box *s1b = (GScroll1Box *) g;

    if ( s1b==NULL )
	return;
    GDrawCancelTimer(s1b->pressed);
    for ( int c=0; c<s1b->count; ++c )
	GGadgetDestroy(s1b->children[c]);
    if ( s1b->sb!=NULL )
	GGadgetDestroy((GGadget *) s1b->sb);
    _ggadget_destroy(g);
}

static void GScroll1BoxMove(GGadget *g, int32 x, int32 y) {
    GScroll1Box *s1b = (GScroll1Box *) g;
    int offx = x-g->r.x, offy = y-g->r.y;

    for ( int c=0; c<s1b->count; ++c ) {
	GGadget *g = s1b->children[c];
	GGadgetMove(g, g->r.x+offx, g->r.y+offy);
    }
    if ( s1b->sb!=NULL )
	GGadgetMove((GGadget *)s1b->sb, s1b->sb->g.r.x+offx, s1b->sb->g.r.y+offy);
    _ggadget_move(g,x,y);
}

/*uint32 count;
    uint32 contentsize;
    uint32 boxsize, boxopposize;
    uint32 minc, curc;
    uint32 ivals;
    uint32 curival;
    uint32 curoffset; */

/* A Scroll1Box will often be embedded in a GHVBox, in which case the desired
 * size will translate in the miminum size of the field. Therefore the desired
 * size should be as small as practical.
 *
 * That means the DesiredSize in that dimension should be the maximum of the
 * DesiredSize of each non-gflowbox child in that dimension and the "squashed
 * size" of each gflowbox child in that dimension. 
 *
 * The scrolling direction should be min scroll size (so that the scrollbar can
 * be operational or the combined height of the children, which ever is
 * smaller.
 */

static void GScroll1BoxGetDesiredSize(GGadget *g, GRect *outer, GRect *inner) {
    GScroll1Box *s1b = (GScroll1Box *) g;
    int bp = GBoxBorderWidth(g->base, g->box);
    int width=0, height=0, c;
    GRect rect;

    for ( c=0; c<s1b->count; ++c ) {
	if (GGadgetIsGFlowBox(s1b->children[c]))
	    _GFlowBoxGetDesiredSize(s1b->children[c], NULL, &rect, true);
	else
	    GGadgetGetDesiredSize(s1b->children[c], NULL, &rect);
	printf(" ds rect: %dx%d\n", rect.width, rect.height);
	if ( s1b->horizontal ) {
	    if ( rect.height > height )
		rect.height = height;
	    width += rect.width;
	    if ( c>0 )
		width += s1b->hpad;
	} else {
	    if ( rect.width > width )
		width = rect.width;
	    height += rect.height;
	    if ( c>0 )
		height += s1b->vpad;
	}
    }
    if ( height > s1b->scrollmin )
	height = s1b->scrollmin;

    if ( inner!=NULL ) {
	inner->x = inner->y = 0;
	inner->width = width;
	inner->height = height;
    }
    if ( outer!=NULL ) {
	outer->x = outer->y = 0;
	outer->width = width+2*bp;
	outer->height = height+2*bp;
    }
}

/*
static void GScroll1Box_ScrollTo(GScroll1Box *s1b, int ival) {
    int newoffset;
    if ( s1b->ivals==0 || sib->curival == ival )
	return;
    if ( !GDrawIsVisible(gl->g.base) )
	return;

    int slide = s1b->contentsize - s1b->boxsize;
    if ( ival+1 == s1b->ivals ) {
	newoffset = s1b->contentsize - s1b->boxsize;
    } else {
	newoffset = s1b->curc * ival;
    }

    GDrawForceUpdate(gl->g.base);

    int offsetdiff = newoffset-s1b->curoffset;
    s1b->curoffset = newoffset;
    s1b->curival = ival;
    for ( int c=0; c<s1b->count; ++c ) {
	GGadget *g = s1b->children[c];
	if ( s1b->horizontal )
	    GGadgetMove(g, g->r.x+offsetdiff, g->r.y);
	else
	    GGadgetMove(g, g->r.x, g->r.y+offsetdiff);
    }
    if ( ydiff>=gl->g.inner.height || -ydiff >= gl->g.inner.height )
	_ggadget_redraw(&gl->g);
    else if ( ydiff!=0 || xoff!=0 )
	GDrawScroll(gl->g.base,&gl->g.inner,xoff,ydiff);
    if ( loff!=0 && gl->vsb!=NULL )
	GScrollBarSetPos(&gl->vsb->g,gl->loff);
}*/

static int GScroll1BoxExpose(GWindow pixmap, GGadget *g, GEvent *event) {
    if ( g->state!=gs_invisible )
	GBoxDrawBorder(pixmap,&g->r,g->box,g->state,false);
    return true;
}


static void GScroll1BoxRedraw(GGadget *g) {
    GScroll1Box *s1b = (GScroll1Box *) g;
    if ( s1b->sb!=NULL )
	GGadgetRedraw((GGadget *) (s1b->sb));
    _ggadget_redraw(g);
}

static void GScroll1BoxResize(GGadget *g, int32 width, int32 height ) {
    GScroll1Box *s1b = (GScroll1Box *) g;
    int bp = GBoxBorderWidth(g->base, g->box);
    int c;
    int old_enabled = GDrawEnableExposeRequests(g->base,false);
    int offset=s1b->curoffset;
    GRect rect;

    width -= 2*bp;
    height -= 2*bp;

    g->inner.x = g->r.x + bp;
    g->inner.y = g->r.y + bp;

    if (s1b->horizontal) {
	if ( g->inner.height!=height ) {
	    for ( c=0; c<s1b->count; ++c ) {
		rect.width = -1;
		rect.height = height;
		// Only likely to work well with a GFlowBox
		GGadgetSetDesiredSize(s1b->children[c], NULL, &rect);
		GGadgetGetDesiredSize(s1b->children[c], NULL, &rect);
		GGadgetResize(s1b->children[c], rect.width, height);
		if ( c>0 )
		    offset += s1b->hpad;
		GGadgetMove(s1b->children[c], g->inner.x+offset, g->inner.y);
		offset += rect.width;
	    }
	}
    } else {
	if ( g->inner.width!=width ) {
	    for ( c=0; c<s1b->count; ++c ) {
		rect.height = -1;
		rect.width = width;
		// Only likely to work well with a GFlowBox
		GGadgetSetDesiredSize(s1b->children[c], NULL, &rect);
		GGadgetGetDesiredSize(s1b->children[c], NULL, &rect);
		GGadgetResize(s1b->children[c], width, rect.height);
		if ( c>0 )
		    offset += s1b->vpad;
		GGadgetMove(s1b->children[c], g->inner.x, g->inner.y+offset);
		offset += rect.height;
	    }
	}
    }
    
    _ggadget_resize(g,width,height);
    GDrawEnableExposeRequests(g->base,old_enabled);
    GDrawRequestExpose(g->base,&g->r,false);
}

static void GScroll1BoxSetVisible(GGadget *g, int visible) {
    GScroll1Box *s1b = (GScroll1Box *) g;
    if ( s1b->sb!=NULL )
	GGadgetSetVisible(&s1b->sb->g,visible);
    _ggadget_setvisible(g,visible);
}

static int GScroll1BoxFillsWindow(GGadget *g) {
    return true;
}

struct gfuncs gscroll1box_funcs = {
    0,
    sizeof(struct gfuncs),

    GScroll1BoxExpose,
    _ggadget_noop,
    _ggadget_noop,
    NULL,
    NULL,
    NULL,
    NULL,

    GScroll1BoxRedraw,
    GScroll1BoxMove,
    GScroll1BoxResize,
    GScroll1BoxSetVisible,
    _ggadget_setenabled,
    _ggadget_getsize,
    _ggadget_getinnersize,

    GScroll1BoxDestroy,

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

    GScroll1BoxGetDesiredSize,
    _ggadget_setDesiredSize,
    GScroll1BoxFillsWindow,
    NULL
};

void GScroll1BoxSetPadding(GGadget *g, int hpad, int vpad) {
    GScroll1Box *s1b = (GScroll1Box *) g;
    if ( hpad>=0 )
        s1b->hpad = hpad;
    if ( vpad>=0 )
        s1b->vpad = vpad;
}

GGadget *GScroll1BoxCreate(struct gwindow *base, GGadgetData *gd,void *data) {
    GScroll1Box *s1b = calloc(1, sizeof(GScroll1Box));

    if ( !gscroll1box_inited )
	GScroll1Box_Init();

    for ( s1b->count=0; gd->u.boxelements[s1b->count]!=NULL; ++s1b->count );

    s1b->g.funcs = &gscroll1box_funcs;
    _GGadget_Create(&s1b->g,base,gd,data,&scroll1box_box);
    s1b->hpad = s1b->vpad = GDrawPointsToPixels(base, 2);
    s1b->scrollmin = GDrawPointsToPixels(base, 30);
    //s1b->curoffset = GDrawPointsToPixels(base, 40);

    s1b->g.takes_input = false;
    s1b->g.takes_keyboard = false;
    s1b->g.focusable = false;

    s1b->children = malloc(s1b->count*sizeof(GGadget *));
    for ( int c=0; c<s1b->count; ++c ) {
        GGadgetCreateData *gcd = gd->u.boxelements[c];
        if ( gcd==GCD_HPad10 )
            s1b->children[c] = (GGadget *) gcd;
        else {
            gcd->gd.pos.x = gcd->gd.pos.y = 1;
            s1b->children[c] = gcd->ret = (gcd->creator)(base,&gcd->gd,gcd->data);
            gcd->ret->contained = true;
        }
    }
    return &s1b->g;
}

void GScroll1BoxFitWindow(GGadget *g) {
    GRect outer, cur, screen;

    GScroll1BoxGetDesiredSize(g, &outer, NULL);
    GDrawGetSize(GDrawGetRoot(NULL),&screen);
    if ( outer.width > screen.width-20 ) outer.width = screen.width-20;
    if ( outer.height > screen.height-40 ) outer.height = screen.height-40;
    GDrawGetSize(g->base,&cur);
    /* Make any offset simmetrical */
    outer.width += 2*g->r.x;
    outer.height += 2*g->r.y;
    printf("  cur: %dx%d  outer: %dx%d\n", cur.width, cur.height, outer.width, outer.height);
    if ( cur.width!=outer.width || cur.height!=outer.height ) {
	GDrawResize(g->base, outer.width, outer.height );
        /* We want to get the resize before we set the window visible */
        /*  and window managers make synchronizing an issue... */
	GDrawSync(GDrawGetDisplayOfWindow(g->base));
	GDrawProcessPendingEvents(GDrawGetDisplayOfWindow(g->base));
	GDrawSync(GDrawGetDisplayOfWindow(g->base));
	GDrawProcessPendingEvents(GDrawGetDisplayOfWindow(g->base));
    } else
	GGadgetResize(g, outer.width-2*g->r.x, outer.height-2*g->r.y );
}

GResInfo *_GScroll1BoxRIHead(void) {
    GScroll1Box_Init();
    return &gscroll1box_ri;
}

