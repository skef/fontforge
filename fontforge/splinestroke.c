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
    // For circle pens
    bigreal radius;
    bigreal radius2; // Squared
    enum linejoin join;
    enum linecap cap;
    bigreal miterlimit;	// For miter joins
	    // PostScript uses 1/sin( theta/2 ) as their miterlimit
	    //  I use -cos(theta). (where theta is the angle between the slopes
	    //  same idea, different implementation
    // For polygon pens
    int n;
    BasePoint *corners;	// Expressed as a vector from the stroke point
    BasePoint *utanvecs; // Normalized slope of line approaching point
    unsigned int remove_inner: 1;
    unsigned int remove_outer: 1;
    unsigned int leave_users_center: 1;	// Don't move the pen so its center is at the origin
    unsigned int scaled_or_rotated: 1;
    unsigned int transform_needed: 1;
    unsigned int open: 1;	// Per-contour: whether the contour being processed is open
    real transform[6];
    real inverse[6];
} StrokeContext;

static char *glyphname=NULL;

static void PrintCorners(StrokeContext *c) {
    for (int i = 0; i < c->n; i++)
        printf("Corner %d: x=%lf, y=%lf, ut=%lf,%lf\n", i, c->corners[i].x, c->corners[i].y, c->utanvecs[i].x, c->utanvecs[i].y);
}

static void OrderCorners(StrokeContext *c, int first) {
    int n = c->n;

    if (first == 0) {
	PrintCorners(c);
	return;
    }

    BasePoint *tbp = malloc(n*sizeof(BasePoint));
    BasePoint *tc = malloc(n*sizeof(BasePoint));

    memcpy(tbp, c->corners, n*sizeof(BasePoint));
    memcpy(c->corners, tbp+first, (n-first)*sizeof(BasePoint));
    memcpy(c->corners+(n-first), tbp, first*sizeof(BasePoint));
    memcpy(tc, c->utanvecs, n*sizeof(BasePoint));
    memcpy(c->utanvecs, tc+first, (n-first)*sizeof(BasePoint));
    memcpy(c->utanvecs+(n-first), tc, first*sizeof(BasePoint));
    PrintCorners(c);
    free(tbp);
    free(tc);
}

static int UTanAngleGreater(BasePoint uta, BasePoint utb) {
    if ( UTanVecNear(uta, utb) )
	return 0;

    if (uta.y >= 0) {
	if (utb.y < 0)
	    return 1;
	return uta.x < utb.x;
    }
    if (utb.y >= 0)
	return 0;
    return uta.x > utb.x;
}

/* Returns the index into c->corners and c->utanvecs of the
 * corner with the largest incoming slope less than or 
 * equal to s. If there is no such entry -1 is returned
 */

static int CornerIndex(StrokeContext *c, BasePoint ut, int is_right, int cw) {
    int i;

    // We trace the right paths backwards so reverse the angle to look for
    if (is_right) {
	ut.x *= -1;
	ut.y *= -1;
    }

    for (i = 0; i<c->n; ++i) {
	if ( UTanVecNear(ut, c->utanvecs[i]) )
	    // Clockwise means decreasing angle and therefore i, ccw means increasing.
	    // The goal is to return the corner continuous with the following curve.
	    // So if we're near the angle and the angles are decreasing return the "lower"
	    // corner, otherwise return the higher one.
	    return cw ? i : (c->n+i-1)%c->n;
	if ( UTanAngleGreater(ut, c->utanvecs[i]) )
	    return (c->n+i-1)%c->n;
    }

    return c->n-1;
}

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
/*  4. It must be convex */
/*  5. No extranious points on the edges */

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

/* Copies the portion of orig from t_start to t_end and translates
 * it so that it starts at point start. The copy is then appended after 
 * point start and the new end-point is returned. When t_start is larger
 * than t_end the spline direction will be reversed and the point
 * corresponding to t_end on the original spline will correspond to
 * start on the new one.
 *
 * Calculations cribbed from https://stackoverflow.com/a/879213
 */
SplinePoint *AppendCubicSplinePortion(Spline *s, bigreal t_start, bigreal t_end,
                                 SplinePoint *start) {
    extended u_start = 1-t_start, u_end = 1-t_end;
    SplinePoint *end;
    BasePoint v, qf, qcf, qct, qt;

    // Intermediate calculations
    qf.x = s->from->me.x*u_start*u_start + s->from->nextcp.x*2*t_start*u_start + s->to->prevcp.x*t_start*t_start;
    qcf.x = s->from->me.x*u_end*u_end + s->from->nextcp.x*2*t_end*u_end + s->to->prevcp.x*t_end*t_end;
    qct.x = s->from->nextcp.x*u_start*u_start + s->to->prevcp.x*2*t_start*u_start+s->to->me.x*t_start*t_start;
    qt.x = s->from->nextcp.x*u_end*u_end + s->to->prevcp.x*2*t_end*u_end + s->to->me.x*t_end*t_end;

    qf.y = s->from->me.y*u_start*u_start + s->from->nextcp.y*2*t_start*u_start + s->to->prevcp.y*t_start*t_start;
    qcf.y = s->from->me.y*u_end*u_end + s->from->nextcp.y*2*t_end*u_end + s->to->prevcp.y*t_end*t_end;
    qct.y = s->from->nextcp.y*u_start*u_start + s->to->prevcp.y*2*t_start*u_start+s->to->me.y*t_start*t_start;
    qt.y = s->from->nextcp.y*u_end*u_end + s->to->prevcp.y*2*t_end*u_end + s->to->me.y*t_end*t_end;

    // Difference vector to offset other points
    v.x = start->me.x - (qf.x*u_start + qct.x*t_start);
    v.y = start->me.y - (qf.y*u_start + qct.y*t_start);
    //printf("vx = %lf, vy = %lf\n", v.x, v.y);

    end = SplinePointCreate(qcf.x*u_end + qt.x*t_end + v.x, qcf.y*u_end + qt.y*t_end + v.y);

    start->nonextcp = false; end->noprevcp = false;
    start->nextcp.x = qf.x*u_end + qct.x*t_end + v.x;
    start->nextcp.y = qf.y*u_end + qct.y*t_end + v.y;
    end->prevcp.x = qcf.x*u_start + qt.x*t_start + v.x;
    end->prevcp.y = qcf.y*u_start + qt.y*t_start + v.y;

    SplineMake3(start,end);

    if ( SplineIsLinear(start->next)) { // Linearish?
        start->nextcp = start->me;
        end->prevcp = end->me;
        start->nonextcp = end->noprevcp = true;
        SplineRefigure(start->next);
    } 
    return end;
}

static bigreal PolyNextT(StrokeContext *c, Spline *s, int i, bigreal cur_t, BasePoint *cur_ut, int is_right) {
    int si[3], j, cnt = 2;
    bigreal t, min_t = 2;
    BasePoint ut = *cur_ut, min_ut;

    // We trace the right paths backwards so reverse the angle to compare against
    if (is_right) {
	ut.x *= -1;
	ut.y *= -1;
    }

    // Look for angles at higher and lower angle corners.
    si[0] = (i+1)%c->n;
    si[1] = i;
    // If we're on a corner slope we need to look for corners on either side
    // but also the current corner because of inflection points.
    if ( UTanVecNear(c->utanvecs[i], ut) ) {
	si[2] = (c->n+i-1)%c->n;
	++cnt;
    }
    printf("i: %d, si0: %d, si1: %d, si2: %d\n", i, si[0], si[1], si[2]);

    for ( j=0; j<cnt; ++j ) {
	ut = c->utanvecs[si[j]];
	if (is_right) {
	    ut.x *= -1;
	    ut.y *= -1;
	}
	t = SplineSolveForUTanVec(s, ut, cur_t);
	if ( t != -1 && t < min_t ) {
	    min_ut = ut;
	    min_t = t;
	}
    }

    if ( min_t>=0 && min_t<=1 ) {
	*cur_ut = min_ut;
	return min_t;
    } else
	return -1;
}

static int ShapeCapNext(int n, int start_i, int end_i, int ccw) {
    if ( ccw ) {
	if ( start_i < end_i ) {
	    if ( start_i >= 0 )
		return --start_i;
	    else if ( end_i < n-1 )
		return n-1;
	    else
		return -1;
	} else if ( start_i > end_i+1 )
	    return --start_i;
	else
	    return -1;
    } else {
	if ( start_i > end_i ) {
	    if ( start_i < n-1 )
		return ++start_i;
	    else if ( end_i > 0 )
		return 0;
	    else
		return -1;
	} else if ( start_i < end_i-1 )
	    return ++start_i;
	else
	    return -1;
    }
}

static void ShapeCap(StrokeContext *c, SplinePoint *start_p, int start_i, SplinePoint *end_p, int end_i, int ccw) {
    SplinePoint *tmp_p;
    int pi = start_i, i;

    while ( (i = ShapeCapNext(c->n, pi, end_i, ccw)) != -1 ) {
	tmp_p = SplinePointCreate(start_p->me.x-c->corners[pi].x+c->corners[i].x, start_p->me.y-c->corners[pi].y+c->corners[i].y);
	SplineMake3(start_p, tmp_p);
	start_p = tmp_p;
	pi = i;
    }
    SplineMake3(start_p, end_p);
}

static void ShapeJoint(StrokeContext *c, SplinePoint *start_p, int start_i, SplinePoint *end_p, int end_i, int is_neg, int is_right ) {
    if ( !is_neg == !is_right ) {
	SplineMake3(start_p, end_p);
	return;
    }
    ShapeCap(c, start_p, start_i, end_p, end_i, is_right);
}

static SplineSet *SplinesToContoursPoly(SplineSet *ss, StrokeContext *c) {
    Spline *s, *first=NULL;
    int is_right, cci, ncci, linear, cw_start, cw_mid;
    int fli = -1, fri = -1, li, ri;
    bigreal last_t, t, tx, ty;
    BasePoint ut_ini = { 2,0 }, ut_start, ut_mid, ut_end, ut_diff;
    SplineSet *left=NULL, *right=NULL, *cur;
    SplinePoint *sp;
    Spline *lspline;

    if ( !c->remove_outer || c->open )
	left = chunkalloc(sizeof(SplineSet));
    if ( !c->remove_inner || c->open )
	right = chunkalloc(sizeof(SplineSet));

    for ( s=ss->first->next; s!=NULL && s!=first; s=s->to->next ) {
	if ( first==NULL )
	    first = s;

	if ( AdjustedSplineLength(s)==0 )
	    continue;		/* We can safely ignore it because it is of zero length */

	ut_start = SplineUTanVecAt(s, 0.0);
	if ( ut_ini.x==2 )
	    ut_ini = ut_start;
	linear = SplineIsLinearish(s);
	if ( linear ) {
	    ut_end = ut_start;
	    cw_start = 0;
	} else {
	    ut_end = SplineUTanVecAt(s, 1.0);
	    cw_start = SplineTurningCWAt(s, 0.0);
	}

	// Left then right
	for ( is_right=0; is_right<=1; ++is_right ) {

	    cci = CornerIndex(c, ut_start, is_right, cw_start);

	    if ( is_right ) {
		if ( right==NULL )
		    continue;
		cur = right;
		if ( fri==-1 )
		    fri = cci;
	    } else {
		if ( left==NULL )
		    continue;
		cur = left;
		if ( fli==-1 )
		    fli = cci;
	    }

	    // Create or verify initial spline point
	    tx = s->from->me.x+c->corners[cci].x;
	    ty = s->from->me.y+c->corners[cci].y;
	    printf("is_right=%d cci=%d corners[cci].x=%lf corners[cci].y=%lf x=%lf, y=%lf, ut_start=%lf,%lf, cw_start=%d\n", is_right, cci, c->corners[cci].x, c->corners[cci].y, tx, ty, ut_start.x, ut_start.y, cw_start);
	    if ( cur->first==NULL ) {
		cur->first = SplinePointCreate(tx, ty);
		cur->last = cur->first;
	    } else if ( cur->last->me.x != tx || cur->last->me.y != ty ) {
		ut_diff = UTanVecDiff(ut_end, ut_start);
		sp = SplinePointCreate(tx, ty);
		ShapeJoint(c, cur->last, is_right ? ri : li, sp, cci, ut_diff.y >= 0, is_right);
		cur->last = sp;
	    }

	    // Complete the path
	    if ( linear ) {
		sp = SplinePointCreate(s->to->me.x+c->corners[cci].x, s->to->me.y+c->corners[cci].y);
		SplineMake3(cur->last, sp);
		cur->last = sp;
	    } else {
		t = 0.0;
		ut_mid = ut_start;
		ncci = -1;
		while ( t < 1.0 ) {
		    last_t = t;
		    t = PolyNextT(c, s, cci, t, &ut_mid, is_right);
		    if ( t==-1.0 || t<last_t ) // That slope not found -- copy rest of spline
			t = 1.0;
		    printf("cci: %d, ncci: %d, ut_mid=%lf,%lf, t=%lf\n", cci, ncci, ut_mid.x, ut_mid.y, t);
		    sp = AppendCubicSplinePortion(s, last_t, t, cur->last);
		    cur->last = sp;
		    if (t < 1.0) {
			cw_mid = SplineTurningCWAt(s, t);
			ncci = CornerIndex(c, ut_mid, is_right, cw_mid);
			sp = SplinePointCreate(sp->me.x-c->corners[cci].x+c->corners[ncci].x, sp->me.y-c->corners[cci].y+c->corners[ncci].y);
			SplineMake3(cur->last, sp);
			cur->last = sp;
			cci = ncci;
		    }
		}
	    }

	    // Record the corner position
	    if ( is_right )
		ri = cci;
	    else
		li = cci;
	}
    }
    // If one SplineSet is empty both should be, return NULL
    if ( (left && left->first==NULL) || (right && right->first==NULL) ) { // No non-zero-length splines.
	chunkfree(left, sizeof(SplineSet));
	chunkfree(right, sizeof(SplineSet));
	return NULL;
    }
/*    if ( left && right) {
	left->next = right;
	left->first->name = copy("L");
	return left;
    } */
    if ( c->open ) {
	SplineSetReverse(right);
	ShapeCap(c, left->last, li, right->first, ri, 0);
	ShapeCap(c, right->last, fri, left->first, fli, 0);
	left->last = left->first;
	right = NULL;
    } else {
	if ( left!=NULL ) {
	    // XXX Close equality test
	    SplineMake3(left->last, left->first);
	    left->last = left->first;
	}
	if ( right!=NULL ) {
	    if ( left != NULL ) {
		SplineSetReverse(right);
		left->next = right;
	    } else
		left = right;
	    SplineMake3(right->last, right->first);
	    right->last = right->first;
	    right = NULL;
	}
    }
    return left;
}

#define NIPOINTS 20
#define NIDIFF (1.0/NIPOINTS)

static SplineSet *SplinesToContoursRadial(SplineSet *ss, StrokeContext *c) {
    Spline *s, *first;
    bigreal length;
    int i;
    bigreal t;
    SplineSet *head = NULL, *last = NULL, *cur;
    StrokePoint *p, stp[NIPOINTS];
    TPoint l[NIPOINTS], r[NIPOINTS];
    SplinePoint *sp;

    first = NULL;
    for ( s=ss->first->next; s!=NULL && s!=first; s=s->to->next ) {
	if ( first==NULL ) first = s;
	length = AdjustedSplineLength(s);
	if ( length==0 )
	    continue;		/* We can safely ignore it because it is of zero length */
	for ( i=0, t=0; i<NIPOINTS; ++i, t+= NIDIFF ) {
	    p = &stp[i];
	    if ( i==NIPOINTS-1 ) t = 1.0;	/* In case there were rounding errors */
	    // FindSlope(p,s,t,NIDIFF,i!=0);
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
	//sp->pointtype = pt_corner;
	// SetNextCPSlope(sp, stp[0].slope, 1);
	cur->first = cur->last = sp;

	sp = SplinePointCreate(l[NIPOINTS-1].x, l[NIPOINTS-1].y);
	//sp->pointtype = pt_corner;
	// SetPrevCPSlope(sp, stp[NIPOINTS-1].slope, -1);
	ApproximateSplineFromPointsSlopes(cur->last,sp,l+1,NIPOINTS-2,false);
	cur->last = sp;

	sp = SplinePointCreate(r[NIPOINTS-1].x, r[NIPOINTS-1].y);
	SplineMake3(cur->last,sp);
	sp->pointtype = pt_corner;
	// SetNextCPSlope(sp, stp[NIPOINTS-1].slope, -1);
	cur->last = sp;

	sp = SplinePointCreate(r[0].x, r[0].y);
	sp->pointtype = pt_corner;
	// SetPrevCPSlope(sp, stp[0].slope, 1);
	ApproximateSplineFromPointsSlopes(cur->last,sp,r+1,NIPOINTS-2,false);
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
	if ( c->pentype == pt_circle )
	    ret = SplinesToContoursRadial(base, c);
	else
	    ret = SplinesToContoursPoly(base, c);
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
    int n, i, max, has_max_utanangle;
    bigreal d2, maxd2, len, sn, co, factor;
    DBounds b;
    BasePoint center, max_utanangle;
    real trans[6];

    if ( si->stroke_type==si_centerline )
	IError("centerline not handled");

    memset(&c,0,sizeof(c));
    c.join = si->join;
    c.cap  = si->cap;
    c.miterlimit = /* -cos(theta) */ -.98 /* theta=~11 degrees, PS miterlimit=10 */;
    c.radius = si->radius;
    c.radius2 = si->radius*si->radius;
    c.remove_inner = si->removeinternal;
    c.remove_outer = si->removeexternal;
    c.leave_users_center = si->leave_users_center;
    c.scaled_or_rotated = si->factor!=NULL;
    if ( si->stroke_type==si_std ) {
	c.pentype = pt_circle;
	if ( si->minorradius==0 )
	    si->minorradius = si->radius;
	if ( si->minorradius!=si->radius ) {
	    c.transform_needed = true;
	    sn = sin(si->penangle);
	    co = cos(si->penangle);
	    factor = si->radius/si->minorradius;
	    c.transform[0] = c.transform[3] = co;
	    c.transform[1] = -sn;
	    c.transform[2] = sn;
	    c.transform[1] *= factor;
	    c.transform[3] *= factor;
	    c.inverse[0] = c.inverse[3] = co;
	    c.inverse[1] = sn;
	    c.inverse[2] = -sn;
	    c.inverse[3] /= factor;
	    c.inverse[2] /= factor;
	}
	ret = SplineSets_Stroke(ss,&c,order2);
    } else if ( si->stroke_type==si_caligraphic ) {
	// This is more or less si_poly with a different setup routine
	c.pentype = pt_square;
	memset(trans,0,sizeof(trans));
	printf("Radius: %lf, Minor Radius: %lf\n", si->radius, si->minorradius);
	// If penangle is 0 and the two radii are equal this should just
	// be the identity matrix
	sn = sin(si->penangle);
	co = cos(si->penangle);
	factor = si->radius/si->minorradius;
	trans[0] = trans[3] = co;
	trans[1] = sn;
	trans[2] = -sn;
	trans[3] /= factor;
	trans[2] /= factor;
	c.n = 4;
	c.corners = malloc(c.n*sizeof(BasePoint));
	c.utanvecs = malloc(c.n*sizeof(BasePoint));
	c.corners[0].x = -si->radius; c.corners[0].y = si->radius;
	c.corners[1].x = si->radius; c.corners[1].y = si->radius;
	c.corners[2].x = si->radius; c.corners[2].y = -si->radius;
	c.corners[3].x = -si->radius; c.corners[3].y = -si->radius;
	for ( i=0; i<c.n; ++i )
	    BpTransform(c.corners+i, c.corners+i, trans);
	has_max_utanangle = -1;
	max_utanangle.x = -1; max_utanangle.y = -DBL_MIN;
	for ( i=0; i<c.n; ++i ) {
//	    printf("corners[i].y: %lf corners[(c.n+i-1)%c.n].y: %lf, corners[i].x: %lf, corners[(c.n+i-1)%c.n].x: %lf\n",
//	           c.corners[i].y, c.corners[(c.n+i-1)%c.n].y, c.corners[i].x, c.corners[(c.n+i-1)%c.n].x);
	    c.utanvecs[i] = UTanVectorize(c.corners[i].x - c.corners[(c.n+i-1)%c.n].x, c.corners[i].y - c.corners[(c.n+i-1)%c.n].y);
	    if ( UTanAngleGreater(c.utanvecs[i], max_utanangle) ) {
		has_max_utanangle = i;
		max_utanangle = c.utanvecs[i];
	    }
	}
	OrderCorners(&c, has_max_utanangle);
	ret = SplineSets_Stroke(ss,&c,order2);
	free(c.corners);
	free(c.utanvecs);
    } else if ( si->stroke_type==si_poly ) {
	c.pentype = pt_poly;
	first = last = NULL;
	max = 0;
	memset(&center,0,sizeof(center));
	for ( active=si->poly; active!=NULL; active=active->next ) {
	    for ( sp=active->first, n=0; ; ) {
		++n;
		if ( sp->next==NULL )
		    return( NULL ); // Error: must be closed
		sp=sp->next->to;
		if ( sp==active->first )
	    break;
	    }
	    if ( n>max )
		max = n;
	}
	c.corners = malloc(max*sizeof(BasePoint));
	c.utanvecs  = malloc(max*sizeof(BasePoint));
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
	    maxd2 = 0;
	    has_max_utanangle = -1;
	    max_utanangle.x = -1; max_utanangle.y = -DBL_MIN;
	    for ( sp=active->first, n=0; ; ) {
		nsp = sp->next->to;
		c.corners[n] = nsp->me;
		c.utanvecs[n] = UTanVectorize(nsp->me.x - sp->me.x, nsp->me.y - sp->me.y);
		if ( UTanAngleGreater(c.utanvecs[n], max_utanangle) ) {
		    has_max_utanangle = n;
		    max_utanangle = c.utanvecs[n];
		}
		++n;
		sp=nsp;
		if ( sp==active->first )
		    break;
	    }
	    c.n = n;
	    OrderCorners(&c, has_max_utanangle);
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
	free(c.utanvecs);
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
