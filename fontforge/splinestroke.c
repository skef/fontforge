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

#include "splinestroke.h"

#include "baseviews.h"
#include "cvundoes.h"
#include "fontforge.h"
#include "splinefont.h"
#include "splineorder2.h"
#include "splineoverlap.h"
#include "splineutil.h"
#include "splineutil2.h"

#include <math.h>

#define PI      3.1415926535897932

#define SP_RIGHT 0
#define SP_LEFT 1

typedef struct strokepoint {
    Spline *sp;
    bigreal t;
    BasePoint me;		/* This lies on the original curve */
    BasePoint slope;		/* Slope at that point */
    BasePoint edge[2];
    unsigned int needs_point: 1;
    //uint8 et;
} StrokePoint;

enum pentype { pt_circle, pt_square, pt_poly };

typedef struct strokecontext {
    enum pentype pentype;
    bigreal resolution;	/* take samples roughly this many em-units */ /* radius/16.0? */
    bigreal radius;
    bigreal radius2;	/* Squared */
    enum linejoin join;	/* Only for circles */
    enum linecap cap;	/* Only for circles */
    bigreal miterlimit;	/* Only for circles if join==lj_miter */
	    /* PostScript uses 1/sin( theta/2 ) as their miterlimit */
	    /*  I use -cos(theta). (where theta is the angle between the slopes*/
	    /*  same idea, different implementation */
    int n;		/* For polygon pens */
    BasePoint *corners;	/* Expressed as an offset from the center of the poly */ /* (Where center could mean center of bounding box, or average or all verteces) */
    BasePoint *slopes;	/* slope[0] is unitvector corners[1]-corners[0] */
    bigreal longest_edge;
    unsigned int open: 1;	/* Original is an open contour */
    unsigned int remove_inner: 1;
    unsigned int remove_outer: 1;
    /* unsigned int rotate_relative_to_direction: 1; */	/* Rotate the polygon pen so it maintains the same orientation with respect to the contour's slope */
    /* Um, the above is essentially equivalent to a circular pen. No point in duplicating it */
    unsigned int leave_users_center: 1;			/* Don't move the pen so its center is at the origin */
    unsigned int scaled_or_rotated: 1;
    unsigned int transform_needed: 1;
    real transform[6];
    real inverse[6];
} StrokeContext;

/*
static void dumpBasePoint(FILE *f, char *nm, BasePoint *p) {
	if (p == NULL)
		return;
	fprintf(f, "   <%s x='%f' y='%f'/>\n", nm, p->x, p->y);
}

static void dumpStrokePoint(FILE *f, StrokePoint *p) {
	fprintf(f, "  <strokepoint>\n");
	dumpBasePoint(f, "me", &(p->me));
	dumpBasePoint(f, "edge", &(p->edge));
	fprintf(f, "  </strokepoint>\n");
}

static void dumpStrokeContext(StrokeContext *c) {
	FILE *f;
	if (c == NULL)
		return;
	f = fopen("/tmp/strokecontext.xml", "w");
	if (f == NULL)
		return;

	fprintf(f, "<strokecontext>\n");
	fprintf(f, " <resolution>%f</resolution>\n", c->resolution);
	fprintf(f, " <miterlimit>%f</miterlimit>\n", c->miterlimit);
	fprintf(f, " <longestedge>%f</longestedge>\n", c->longest_edge);
	fprintf(f, " <open>%d</open>\n", (int) c->open);
	for (int k = 0; k < c->cur; k++) {
		dumpStrokePoint(f, &(c->all[k]));
	}
	fprintf(f, "</strokecontext>\n");
	fclose(f);
}
*/

static char *glyphname=NULL;

static int AdjustedSplineLength(Spline *s) {
    /* I find the spline length in hopes of subdividing the spline into chunks*/
    /*  roughly 1em unit apart. But if the spline bends a lot, then when      */
    /*  expanded out by a pen the distance between chunks will be amplified   */
    /*  So do something quick and dirty to adjust for that */
    /* What I really want to do is to split the spline into segments, find */
    /*  the radius of curvature, find the length, multiply length by       */
    /*  (radius-curvature+c->radius)/radius-curvature                      */
    /*  For each segment and them sum that. But that's too hard */
    bigreal len = SplineLength(s);
    bigreal xdiff=s->to->me.x-s->from->me.x, ydiff=s->to->me.y-s->from->me.y;
    bigreal distance = sqrt(xdiff*xdiff+ydiff*ydiff);

    if ( len<=distance )
return( len );			/* It's a straight line */
    len += 1.5*(len-distance);
return( len );
}

/* Something is a valid polygonal pen if: */
/*  1. It contains something (not a line, or a single point, something with area) */
/*  2. It contains ONE contour */
/*	(Multiple contours can be dealt with as long as each contour follows */
/*	 requirements 3-n. If we have multiple contours: Find the offset from */
/*	 the center of each contour to the center of all the contours. Trace  */
/*	 each contour individually and then translate by this offset) */
/*  3. All edges must be lines (no control points) */
/*  4. No more than 255 corners (could extend this number, but why?) */
/*  5. It must be drawn clockwise (not really an error, if found just invert, but important for later checks). */
/*  6. It must be convex */
/*  7. No extranious points on the edges */

enum PolyType PolygonIsConvex(BasePoint *poly,int n, int *badpointindex) {
    /* For each vertex: */
    /*  Remove that vertex from the polygon, and then test if the vertex is */
    /*  inside the resultant poly. If it is inside, then the polygon is not */
    /*  convex */
    /* If all verteces are outside then we have a convex poly */
    /* If one vertex is on an edge, then we have a poly with an unneeded vertex */
    bigreal nx,ny;
    int i,j,ni;

    if ( badpointindex!=NULL )
	*badpointindex = -1;
    if ( n<3 )
return( Poly_TooFewPoints );
    /* All the points might lie on one line. That wouldn't be a polygon */
    nx = -(poly[1].y-poly[0].y);
    ny = (poly[1].x-poly[0].x);
    for ( i=2; i<n; ++i ) {
	if ( (poly[i].x-poly[0].x)*ny - (poly[i].y-poly[0].y)*nx != 0 )
    break;
    }
    if ( i==n )
return( Poly_Line );		/* Colinear */
    if ( n==3 ) {
	/* Triangles are always convex */
return( Poly_Convex );
    }

    for ( j=0; j<n; ++j ) {
	/* Test to see if poly[j] is inside the polygon poly[0..j-1,j+1..n] */
	/* Basically the hit test code above modified to ignore poly[j] */
	int outside = 0, zero_cnt=0, sign=0;
	bigreal sx, sy, dx,dy, dot;

	for ( i=0; ; ++i ) {
	    if ( i==j )
	continue;
	    ni = i+1;
	    if ( ni==n ) ni=0;
	    if ( ni==j ) {
		++ni;
		if ( ni==n )
		    ni=0;		/* Can't be j, because it already was */
	    }

	    sx = poly[ni].x - poly[i].x;
	    sy = poly[ni].y - poly[i].y;
	    dx = poly[j ].x - poly[i].x;
	    dy = poly[j ].y - poly[i].y;
	    /* Only care about signs, so I don't need unit vectors */
	    dot = dx*sy - dy*sx;
	    if ( dot==0 )
		++zero_cnt;
	    else if ( sign==0 )
		sign= dot;
	    else if ( (dot<0 && sign>0) || (dot>0 && sign<0)) {
		outside = true;
	break;
	    }
	    if ( ni==0 )
	break;
	}
	if ( !outside ) {
	    if ( badpointindex!=NULL )
		*badpointindex = j;
	    if ( zero_cnt>0 )
return( Poly_PointOnEdge );
	    else
return( Poly_Concave );
	}
    }
return( Poly_Convex );
}

static void FindSlope(StrokePoint *p, Spline *s, bigreal t, bigreal tdiff, int peekback) {
    bigreal len;

    p->sp = s;
    p->t = t;
    p->me.x = ((s->splines[0].a*t+s->splines[0].b)*t+s->splines[0].c)*t+s->splines[0].d;
    p->me.y = ((s->splines[1].a*t+s->splines[1].b)*t+s->splines[1].c)*t+s->splines[1].d;
    p->slope.x = (3*s->splines[0].a*t+2*s->splines[0].b)*t+s->splines[0].c;
    p->slope.y = (3*s->splines[1].a*t+2*s->splines[1].b)*t+s->splines[1].c;

    if ( p->slope.x==0 && p->slope.y==0 ) {
	/* If the control point is at the endpoint then at the endpoints we */
	/*  have an undefined slope. Can't have that. */
	/* I suppose it could happen elsewhere */
	if ( peekback )
	    p->slope = p[-1].slope;
	else {
	    bigreal nextt=t+tdiff;
	    p->slope.x = (3*s->splines[0].a*nextt+2*s->splines[0].b)*nextt+s->splines[0].c;
	    p->slope.y = (3*s->splines[1].a*nextt+2*s->splines[1].b)*nextt+s->splines[1].c;
	    if ( p->slope.x==0 && p->slope.y==0 ) {
		BasePoint next;
		next.x = ((s->splines[0].a*nextt+s->splines[0].b)*nextt+s->splines[0].c)*nextt+s->splines[0].d;
		next.y = ((s->splines[1].a*nextt+s->splines[1].b)*nextt+s->splines[1].c)*nextt+s->splines[1].d;
		p->slope.x = next.x - p->me.x;
		p->slope.y = next.y - p->me.y;
	    }
	}
	if ( p->slope.x==0 && p->slope.y==0 ) {
	    p->slope.x = s->to->me.x = s->from->me.x;
	    p->slope.y = s->to->me.y = s->from->me.y;
	}
	if ( p->slope.x==0 && p->slope.y==0 )
	    p->slope.x = 1;
    }
    len = p->slope.x*p->slope.x + p->slope.y*p->slope.y;
    if ( len!=0 ) {
	len = sqrt(len);
	p->slope.x/=len;
	p->slope.y/=len;
    }
}

#define NIPOINTS 20
#define NIDIFF (1.0/NIPOINTS)

static SplineSet *SplinesToContours(SplineSet *ss, StrokeContext *c) {
    Spline *s, *first;
    bigreal length;
    int i, len;
    bigreal diff, t;
    int open = ss->first->prev == NULL;
    int gothere = false;
    SplineSet *head = NULL, *last = NULL, *cur;
    StrokePoint *p, stp[NIPOINTS];
    TPoint l[NIPOINTS], r[NIPOINTS];
    SplinePoint *sp;
    Spline *lspline;

    first = NULL;
    for ( s=ss->first->next; s!=NULL && s!=first; s=s->to->next ) {
	if ( first==NULL ) first = s;
	length = AdjustedSplineLength(s);
	if ( length==0 )		/* This can happen when we have a spline with the same first and last point and no control point */
    continue;		/* We can safely ignore it because it is of zero length */
	for ( i=0, t=0; i<NIPOINTS; ++i, t+= NIDIFF ) {
	    p = &stp[i];
	    if ( i==NIPOINTS-1 ) t = 1.0;	/* In case there were rounding errors */
	    FindSlope(p,s,t,NIDIFF,i!=0);
	    l[i].x = p->me.x - c->radius*p->slope.y;
	    l[i].y = p->me.y + c->radius*p->slope.x;
	    l[i].t = t;
	    printf(" x: %f, y: %f\n", l[i].x, l[i].y);
	    r[i].x = p->me.x + c->radius*p->slope.y;
	    r[i].y = p->me.y - c->radius*p->slope.x;
	    r[i].t = 1.0-t;
	}
	cur = chunkalloc(sizeof(SplineSet));
	if ( head==NULL )
	    head = cur;
	else
	    last->next = cur;
	last = cur;
	sp = SplinePointCreate(l[0].x, l[0].y);
	sp->nextcp.x += stp[0].slope.x;
	sp->nextcp.y += stp[0].slope.y;
	cur->first = cur->last = sp;
	sp = SplinePointCreate(l[NIPOINTS-1].x, l[NIPOINTS-1].y);
	sp->prevcp.x -= stp[NIPOINTS-1].slope.x;
	sp->prevcp.y -= stp[NIPOINTS-1].slope.y;
	ApproximateSplineFromPoints(cur->last,sp,l+1,NIPOINTS-2,false);
	cur->last = sp;
	sp = SplinePointCreate(r[NIPOINTS-1].x, r[NIPOINTS-1].y);
	SplineMake3(cur->last,sp);
	cur->last = sp;
	sp = SplinePointCreate(r[0].x, r[0].y);
	ApproximateSplineFromPoints(cur->last,sp,r+1,NIPOINTS-2,false);
	cur->last = sp;
	SplineMake3(cur->last,cur->first);
	cur->last = cur->first;
    }
    return head;
}

/******************************************************************************/
/* ******************************* To Splines ******************************* */
/******************************************************************************/

static BasePoint SquareCorners[] = {
    { -1,  1 },
    {  1,  1 },
    {  1, -1 },
    { -1, -1 },
};

static struct shapedescrip {
    BasePoint me, prevcp, nextcp;
}
unitcircle[] = {
    { { -1, 0 }, { -1, -0.552 }, { -1, 0.552 } },
    { { 0 , 1 }, { -0.552, 1 }, { 0.552, 1 } },
    { { 1, 0 }, { 1, 0.552 }, { 1, -0.552 } },
    { { 0, -1 }, { 0.552, -1 }, { -0.552, -1 } },
    { { 0, 0 }, { 0, 0 }, { 0, 0 } }
};

static SplinePoint *SpOnCircle(int i,bigreal radius,BasePoint *center) {
    SplinePoint *sp = SplinePointCreate(unitcircle[i].me.x*radius + center->x,
					unitcircle[i].me.y*radius + center->y);
    sp->pointtype = pt_curve;
    sp->prevcp.x = unitcircle[i].prevcp.x*radius + center->x;
    sp->prevcp.y = unitcircle[i].prevcp.y*radius + center->y;
    sp->nextcp.x = unitcircle[i].nextcp.x*radius + center->x;
    sp->nextcp.y = unitcircle[i].nextcp.y*radius + center->y;
    sp->nonextcp = sp->noprevcp = false;
return( sp );
}

SplineSet *UnitShape(int n) {
    SplineSet *ret;
    SplinePoint *sp1, *sp2;
    int i;
    BasePoint origin;

    ret = chunkalloc(sizeof(SplineSet));
    if ( n>=3 || n<=-3 ) {
	/* Regular n-gon with n sides */
	/* Inscribed in a unit circle, if n<0 then circumscribed around */
	bigreal angle = 2*PI/(2*n);
	bigreal factor=1;
	if ( n<0 ) {
	    angle = -angle;
	    n = -n;
	    factor = 1/cos(angle);
	}
	angle -= PI/2;
	ret->first = sp1 = SplinePointCreate(factor*cos(angle), factor*sin(angle));
	sp1->pointtype = pt_corner;
	for ( i=1; i<n; ++i ) {
	    angle = 2*PI/(2*n) + i*2*PI/n - PI/2;
	    sp2 = SplinePointCreate(factor*cos(angle),factor*sin(angle));
	    sp2->pointtype = pt_corner;
	    SplineMake3(sp1,sp2);
	    sp1 = sp2;
	}
	SplineMake3(sp1,ret->first);
	ret->last = ret->first;
	SplineSetReverse(ret);		/* Drat, just drew it counter-clockwise */
    } else if ( n ) {
	ret->first = sp1 = SplinePointCreate(SquareCorners[0].x,
			    SquareCorners[0].y);
	sp1->pointtype = pt_corner;
	for ( i=1; i<4; ++i ) {
	    sp2 = SplinePointCreate(SquareCorners[i].x,
				SquareCorners[i].y);
	    sp2->pointtype = pt_corner;
	    SplineMake3(sp1,sp2);
	    sp1 = sp2;
	}
	SplineMake3(sp1,ret->first);
	ret->last = ret->first;
    } else {
	/* Turn into a circle */
	origin.x = origin.y = 0;
	ret->first = sp1 = SpOnCircle(0,1,&origin);
	for ( i=1; i<4; ++i ) {
	    sp2 = SpOnCircle(i,1,&origin);
	    SplineMake3(sp1,sp2);
	    sp1 = sp2;
	}
	SplineMake3(sp1,ret->first);
	ret->last = ret->first;
    }
return( ret );
}

static SplinePointList *SinglePointStroke(SplinePoint *sp,struct strokecontext *c) {
    SplineSet *ret;
    SplinePoint *sp1, *sp2;
    int i;

    ret = chunkalloc(sizeof(SplineSet));

    if ( c->pentype==pt_circle && c->cap==lc_butt ) {
	/* Leave as a single point */
	ret->first = ret->last = SplinePointCreate(sp->me.x,sp->me.y);
	ret->first->pointtype = pt_corner;
    } else if ( c->pentype==pt_circle && c->cap==lc_round ) {
	/* Turn into a circle */
	ret->first = sp1 = SpOnCircle(0,c->radius,&sp->me);
	for ( i=1; i<4; ++i ) {
	    sp2 = SpOnCircle(i,c->radius,&sp->me);
	    SplineMake3(sp1,sp2);
	    sp1 = sp2;
	}
	SplineMake3(sp1,ret->first);
	ret->last = ret->first;
    } else if ( c->pentype==pt_circle || c->pentype==pt_square ) {
	ret->first = sp1 = SplinePointCreate(sp->me.x+c->radius*SquareCorners[0].x,
			    sp->me.y+c->radius*SquareCorners[0].y);
	sp1->pointtype = pt_corner;
	for ( i=1; i<4; ++i ) {
	    sp2 = SplinePointCreate(sp->me.x+c->radius*SquareCorners[i].x,
				sp->me.y+c->radius*SquareCorners[i].y);
	    sp2->pointtype = pt_corner;
	    SplineMake3(sp1,sp2);
	    sp1 = sp2;
	}
	SplineMake3(sp1,ret->first);
	ret->last = ret->first;
    } else {
	ret->first = sp1 = SplinePointCreate(sp->me.x+c->corners[0].x,
			    sp->me.y+c->corners[0].y);
	sp1->pointtype = pt_corner;
	for ( i=1; i<c->n; ++i ) {
	    sp2 = SplinePointCreate(sp->me.x+c->corners[i].x,
				sp->me.y+c->corners[i].y);
	    sp2->pointtype = pt_corner;
	    SplineMake3(sp1,sp2);
	    sp1 = sp2;
	}
	SplineMake3(sp1,ret->first);
	ret->last = ret->first;
    }
    return( ret );
}

static SplineSet *SplineSet_Stroke(SplineSet *ss,struct strokecontext *c,
	int order2) {
    SplineSet *base;
    SplineSet *ret;

    base = SplinePointListCopy(ss);
    base = SSRemoveZeroLengthSplines(base); /* Hard to get a slope for a zero length spline, that causes problems */
    if ( base==NULL )
	return(NULL);
    if ( c->transform_needed )
	base = SplinePointListTransform(base,c->transform,tpt_AllPoints);
    if ( base->first->next==NULL )
	ret = SinglePointStroke(base->first,c);
    else {
	c->open = base->first->prev==NULL;
	//dumpStrokeContext(c);
	ret = SplinesToContours(base, c);
    }
    if ( c->transform_needed )
	ret = SplinePointListTransform(ret,c->inverse,tpt_AllPoints);
    if ( order2 )
	ret = SplineSetsConvertOrder(ret,order2 );
    SplinePointListFree(base);
    return(ret);
}

static SplineSet *SplineSets_Stroke(SplineSet *ss,struct strokecontext *c,
	int order2) {
    SplineSet *first=NULL, *last=NULL, *cur;

    while ( ss!=NULL ) {
	cur = SplineSet_Stroke(ss,c,order2);
	if ( first==NULL )
	    first = last = cur;
	else
	    last->next = cur;
	if ( last!=NULL )
	    while ( last->next!=NULL )
		last = last->next;
	ss = ss->next;
    }
return( first );
}

SplineSet *SplineSetStroke(SplineSet *ss,StrokeInfo *si, int order2) {
    StrokeContext c;
    SplineSet *first, *last, *cur, *ret, *active = NULL, *anext;
    SplinePoint *sp, *nsp;
    int n, max;
    bigreal d2, maxd2, len, maxlen;
    DBounds b;
    BasePoint center;
    real trans[6];

    if ( si->stroke_type==si_centerline )
	IError("centerline not handled");

    memset(&c,0,sizeof(c));
    c.resolution = si->resolution;
    if ( si->resolution==0 )
	c.resolution = 1;
    c.pentype = si->stroke_type==si_std ? pt_circle : pt_poly;
    c.join = si->join;
    c.cap  = si->cap;
    c.miterlimit = /* -cos(theta) */ -.98 /* theta=~11 degrees, PS miterlimit=10 */;
    c.radius = si->radius;
    c.radius2 = si->radius*si->radius;
    c.remove_inner = si->removeinternal;
    c.remove_outer = si->removeexternal;
    c.leave_users_center = si->leave_users_center;
    c.scaled_or_rotated = si->factor!=NULL;
    if ( c.pentype==pt_circle || c.pentype==pt_square ) {
	if ( si->minorradius==0 )
	    si->minorradius = si->radius;
	if ( si->minorradius!=si->radius ||
		(si->penangle!=0 && si->stroke_type!=si_std) ) {	/* rotating a circle is irrelevant (rotating an elipse means something) */
	    bigreal sn,co,factor;
	    c.transform_needed = true;
	    sn = sin(si->penangle);
	    co = cos(si->penangle);
	    factor = si->radius/si->minorradius;
	    c.transform[0] = c.transform[3] = co;
	    c.transform[1] = -sn;
	    c.transform[2] = sn;
	    c.transform[1] *= factor; c.transform[3] *= factor;
	    c.inverse[0] = c.inverse[3] = co;
	    c.inverse[1] = sn;
	    c.inverse[2] = -sn;
	    c.inverse[3] /= factor; c.inverse[2] /= factor;
	}
	if ( si->resolution==0 && c.resolution>c.radius/3 )
	    c.resolution = c.radius/3;
    if (c.resolution == 0) {
        ff_post_notice(_("Invalid stroke parameters"), _("Stroke resolution is zero"));
        return SplinePointListCopy(ss);
    }
	ret = SplineSets_Stroke(ss,&c,order2);
    } else {
	first = last = NULL;
	max = 0;
	memset(&center,0,sizeof(center));
        // if ( c.pentype==pt_square ) { # make poly for square
	for ( active=si->poly; active!=NULL; active=active->next ) {
	    for ( sp=active->first, n=0; ; ) {
		++n;
		if ( sp->next==NULL )
return( NULL );				/* That's an error, must be closed */
		sp=sp->next->to;
		if ( sp==active->first )
	    break;
	    }
	    if ( n>max )
		max = n;
	}
	c.corners = malloc(max*sizeof(BasePoint));
	c.slopes  = malloc(max*sizeof(BasePoint));
	memset(trans,0,sizeof(trans));
	trans[0] = trans[3] = 1;
	if ( !c.leave_users_center ) {
	    SplineSetQuickBounds(si->poly,&b);
	    trans[4] = -(b.minx+b.maxx)/2;
	    trans[5] = -(b.miny+b.maxy)/2;
	    SplinePointListTransform(si->poly,trans,tpt_AllPoints);
	}
	for ( active=si->poly; active!=NULL; active=anext ) {
	    int reversed = false;
	    if ( !SplinePointListIsClockwise(active) ) {
		reversed = true;
		SplineSetReverse(active);
	    }
	    if ( !c.scaled_or_rotated ) {
		anext = active->next; active->next = NULL;
		SplineSetQuickBounds(active,&b);
		trans[4] = -(b.minx+b.maxx)/2;
		trans[5] = -(b.miny+b.maxy)/2;
		SplinePointListTransform(active,trans,tpt_AllPoints);	/* Only works if pen is fixed and does not rotate or get scaled */
		active->next = anext;
	    }
	    maxd2 = 0; maxlen = 0;
	    for ( sp=active->first, n=0; ; ) {
		nsp = sp->next->to;
		c.corners[n] = sp->me;
		c.slopes[n].x = nsp->me.x - sp->me.x;
		c.slopes[n].y = nsp->me.y - sp->me.y;
		len = c.slopes[n].x*c.slopes[n].x + c.slopes[n].y*c.slopes[n].y;
		len = sqrt(len);
		if ( len>maxlen ) maxlen = len;
		if ( len!=0 ) {
		    c.slopes[n].x/=len; c.slopes[n].y/=len;
		}
		d2 = sp->me.x*sp->me.x + sp->me.y*sp->me.y;
		if ( d2>maxd2 )
		    maxd2 = d2;
		++n;
		sp=nsp;
		if ( sp==active->first )
	    break;
	    }
	    c.n = n;
	    c.longest_edge = maxlen;
	    c.radius = sqrt(maxd2);
	    c.radius2 = maxd2;
	    if ( si->resolution==0 && c.resolution>c.radius/3 )
		c.resolution = c.radius/3;
        if (c.resolution == 0) {
            ff_post_notice(_("Invalid stroke parameters"), _("Stroke resolution is zero"));
            return SplinePointListCopy(ss);
        }
	    cur = SplineSets_Stroke(ss,&c,order2);
	    if ( !c.scaled_or_rotated ) {
		trans[4] = -trans[4]; trans[5] = -trans[5];
		SplinePointListTransform(cur,trans,tpt_AllPoints);
		anext = active->next; active->next = NULL;
		SplinePointListTransform(active,trans,tpt_AllPoints);
		active->next = anext;
	    }
	    if ( reversed ) {
		SplineSet *ss;
		for ( ss=cur; ss!=NULL; ss=ss->next )
		    SplineSetReverse(ss);
		SplineSetReverse(active);
	    }
	    if ( first==NULL )
		first = last = cur;
	    else
		last->next = cur;
	    if ( last!=NULL )
		while ( last->next!=NULL )
		    last = last->next;
	}
	free(c.corners);
	free(c.slopes);
	ret = first;
    }
return( ret );
}

void FVStrokeItScript(void *_fv, StrokeInfo *si,int pointless_argument) {
    FontViewBase *fv = _fv;
    int layer = fv->active_layer;
    SplineSet *temp;
    int i, cnt=0, gid;
    SplineChar *sc;

    for ( i=0; i<fv->map->enccount; ++i ) if ( (gid=fv->map->map[i])!=-1 && fv->sf->glyphs[gid]!=NULL && fv->selected[i] )
	++cnt;
    ff_progress_start_indicator(10,_("Stroking..."),_("Stroking..."),0,cnt,1);

    SFUntickAll(fv->sf);
    for ( i=0; i<fv->map->enccount; ++i ) {
	if ( (gid=fv->map->map[i])!=-1 && (sc = fv->sf->glyphs[gid])!=NULL &&
		!sc->ticked && fv->selected[i] ) {
	    sc->ticked = true;
	    glyphname = sc->name;
	    if ( sc->parent->multilayer ) {
		SCPreserveState(sc,false);
		for ( layer = ly_fore; layer<sc->layer_cnt; ++layer ) {
		    temp = SplineSetStroke(sc->layers[layer].splines,si,sc->layers[layer].order2);
		    SplinePointListsFree( sc->layers[layer].splines );
		    sc->layers[layer].splines = temp;
		}
		SCCharChangedUpdate(sc,ly_all);
	    } else {
		SCPreserveLayer(sc,layer,false);
		temp = SplineSetStroke(sc->layers[layer].splines,si,sc->layers[layer].order2);
		SplinePointListsFree( sc->layers[layer].splines );
		sc->layers[layer].splines = temp;
		SCCharChangedUpdate(sc,layer);
	    }
	    if ( !ff_progress_next())
    break;
	}
    }
 glyphname = NULL;
    ff_progress_end_indicator();
}
