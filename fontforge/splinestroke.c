/* Copyright (C) 2000-2012 by George Williams, 2019 by Skef Iterum */
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
#include <assert.h>

#define PI      3.1415926535897932

#define BASEPOINT_LENSQ(v) ((v).x*(v).x+(v).y*(v).y)
#define BP_REV(v) (BasePoint) { -(v).x, -(v).y }
#define BP_REV_IF(t, v) (t ? (BasePoint) { -(v).x, -(v).y } : (v))
#define BP_ADD(v1, v2) (BasePoint) { (v1).x + (v2).x, (v1).y + (v2).y } 

enum pentype { pt_circle, pt_square, pt_poly };

/* c->nibcorners is a structure that models the "corners" of the current nib
 * treated as if all splines were linear, and therefore as a convex polygon.
 */

typedef struct nibcorner {
    SplinePoint *on_nib;
    BasePoint utanvec;    // Unit tangent of the "line" from the previous point
    unsigned int linear: 1;
} NibCorner;

typedef struct strokecontext {
    enum pentype pentype;
    // For circle pens
    bigreal radius;
    enum linejoin join;
    enum linecap cap;
    bigreal miterlimit;	// For miter joins
	    // PostScript uses 1/sin( theta/2 ) as their miterlimit
	    //  I use -cos(theta). (where theta is the angle between the slopes
	    //  same idea, different implementation
    SplineSet *nib;
    int n;
    NibCorner *nibcorners;
    unsigned int remove_inner: 1;
    unsigned int remove_outer: 1;
    unsigned int leave_users_center: 1;	// Don't move the pen so its center is at the origin
    unsigned int scaled_or_rotated: 1;
} StrokeContext;

#define NIBOFF_HIGH_IDX 0
#define NIBOFF_LOW_IDX 1

// The corner of the nib corresponding to utanvec and the spline
// offset corresponding to that angle
typedef struct niboffset {
    BasePoint utanvec;
    int nci[2];
    BasePoint off[2];
    bigreal nt[2];
    unsigned int linear: 1;  // The nib edge is a line
    unsigned int at_line: 1; // utanvec equals the corner angle (implies linear)
} NibOffset;

static char *glyphname=NULL;

static BasePoint UTanVecDiff(BasePoint ut_ref, BasePoint ut_vec) {
    BasePoint ret;
    ret.x = ut_ref.x*ut_vec.x + ut_ref.y*ut_vec.y; // dot product
    ret.y = ut_ref.x*ut_vec.y - ut_ref.y*ut_vec.x; // mag of cross product
    return ret;
}

static int UTanAngleGreater(BasePoint uta, BasePoint utb) {
    assert(RealNear(BASEPOINT_LENSQ(uta),1) && RealNear(BASEPOINT_LENSQ(utb),1));

    if ( BPNEAR(uta, utb) )
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

static NibOffset *CalcNibOffset(StrokeContext *c, BasePoint ut, NibOffset *no, int i_hint) {
    int i, pi;
    Spline *ns;

    if ( no==NULL )
	no = malloc(sizeof(NibOffset));

    memset(no,0,sizeof(NibOffset));
    no->utanvec = ut;

    // Easier to let the Near case fall back to the search
    if (   i_hint != -1
        && UTanAngleGreater(c->nibcorners[i_hint].utanvec, ut)
        && (   (i_hint+1==c->n)
            || UTanAngleGreater(ut, c->nibcorners[i_hint+1].utanvec)))
	i = i_hint; // Hint was accurate
    else {
	for (i = 0; i<c->n; ++i) {
	    if ( BPNEAR(no->utanvec, c->nibcorners[i].utanvec) ) {
		no->at_line = c->nibcorners[i].linear;
		break;
	    } else if ( UTanAngleGreater(no->utanvec, c->nibcorners[i].utanvec) ) {
		// Go back one
		i = (c->n+i-1)%c->n;
		break;
	    }
	}
    }

    if ( i==c->n )
	--i;

    // i is now the index to corner with least utanvec >= ut (or c->n-1 if
    // ut is larger or smaller than all corners)

    if ( c->nibcorners[i].linear ) {
	no->linear = true;
	no->off[NIBOFF_HIGH_IDX] = c->nibcorners[i].on_nib->me;
	no->nt[NIBOFF_HIGH_IDX] = 1.0;
	if ( no->at_line ) {
	    no->nci[NIBOFF_LOW_IDX] = (c->n+i-1)%c->n;
	    no->off[NIBOFF_LOW_IDX] = c->nibcorners[no->nci[NIBOFF_LOW_IDX]].on_nib->me;
	    no->nt[NIBOFF_LOW_IDX] = 0.0;
	} else
	    no->off[NIBOFF_LOW_IDX] = no->off[NIBOFF_HIGH_IDX];
	    no->nt[NIBOFF_LOW_IDX] = no->nt[NIBOFF_HIGH_IDX];
    } else {
	ns = c->nibcorners[i].on_nib->prev;
	// Nib splines are locally convex and should therefore have t value per slope
	if ( BPNEAR(no->utanvec, c->nibcorners[i].utanvec) )
	    no->nt[0] = no->nt[1] = 0.0;
	else
	    no->nt[0] = no->nt[1] = SplineSolveForUTanVec(ns, no->utanvec, 0.0);
	no->off[0] = no->off[1] = SPLINEPVAL(ns, no->nt[NIBOFF_HIGH_IDX]);
	printf("nt=%lf, off.x=%lf, off.y=%lf\n", no->nt[NIBOFF_HIGH_IDX], no->off[NIBOFF_HIGH_IDX].x, no->off[NIBOFF_HIGH_IDX].y);
    }
    return no;
}

/* A contour is a valid convex polygon if:
 *  1. It contains something (not a line or a single point
 *     but a closed piecewise spline containing area) 
 *  2. All edges/splines are lines
 *  3. It is convex
 *  4. There are no extranious points on the edges
 *
 * A contour is a valid nib shape if:
 *  5. It's on-curve points would form a valid convex
 *     polygon if there were no control points.
 *  6. Every edge/spline is either:
 *     6a. A line, or
 *     6b. A spline where *both* control points are
 *         outside of the polygon but inside the 
 *         triangle formed by the spline edge and
 *         extending the lines of the adjacent edges.
 *
 * A valid nib is a SplineSet containing only valid
 * nib shapes.
 */

enum ShapeType NibIsValid(SplineSet *ss) {
    /* For each vertex: */
    /*  Remove that vertex from the polygon, and then test if the vertex is */
    /*  inside the resultant poly. If it is inside, then the polygon is not */
    /*  convex */
    /* If all verteces are outside then we have a convex poly */
    /* If one vertex is on an edge, then we have a poly with an unneeded vertex */
/*
    bigreal nx,ny;
    int i,j,ni;

    if ( n<3 )
return( Shape_TooFewPoints );
    // All the points might lie on one line. That wouldn't be a polygon
    nx = -(poly[1].y-poly[0].y);
    ny = (poly[1].x-poly[0].x);
    for ( i=2; i<n; ++i ) {
	if ( (poly[i].x-poly[0].x)*ny - (poly[i].y-poly[0].y)*nx != 0 )
    break;
    }
    if ( i==n )
return( Shape_Line );		// Colinear
    if ( n==3 ) {
	// Triangles are always convex
return( Shape_Convex );
    }

    for ( j=0; j<n; ++j ) {
	// Test to see if poly[j] is inside the polygon poly[0..j-1,j+1..n]
	// Basically the hit test code above modified to ignore poly[j]
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
		    ni=0;		// Can't be j, because it already was
	    }

	    sx = poly[ni].x - poly[i].x;
	    sy = poly[ni].y - poly[i].y;
	    dx = poly[j ].x - poly[i].x;
	    dy = poly[j ].y - poly[i].y;
	    // Only care about signs, so I don't need unit vectors
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
//	    if ( badpointindex!=NULL )
//		*badpointindex = j;
	    if ( zero_cnt>0 )
return( Shape_PointOnEdge );
	    else
return( Shape_Concave );
	}
    }
    */
return( Shape_Convex );
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

    if ( RealWithin(t_start, t_end, 1e-3) )
	return start;

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

    if ( SplineIsLinear(start->next)) { // Linearish instead?
        start->nextcp = start->me;
        end->prevcp = end->me;
        start->nonextcp = end->noprevcp = true;
        SplineRefigure(start->next);
    } 
    return end;
}

// src_start and src_end must be on same contour
SplinePoint *AppendCubicSplineSetPortion(SplinePoint *src_start, bigreal t_start,
                                         SplinePoint *src_end, bigreal t_end,
					 SplinePoint *dst_start, int backward) {
    Spline *s = backward ? src_start->prev : src_start->next;

    // Handle single splice case
    if ( src_start==src_end && (( t_start<t_end && !backward ) || (t_start>t_end && backward))) {
	dst_start = AppendCubicSplinePortion(s, t_start, t_end, dst_start);
	return dst_start;
    }

    dst_start = AppendCubicSplinePortion(s, t_start, backward ? 0 : 1, dst_start);

    // walk spline-wise until pointing to src_end
    while ( (backward ? s->from : s->to) != src_end ) {
	assert( (backward ? s->from : s->to)!=src_start );
        s = backward ? s->from->prev : s->to->next;
	assert( s!=NULL );
        dst_start = AppendCubicSplinePortion(s, backward ? 1 : 0, backward ? 0 : 1, dst_start);
    }
    s = backward ? s->from->prev : s->to->next;
    dst_start = AppendCubicSplinePortion(s, backward ? 1 : 0, t_end, dst_start);
    return dst_start;
}

#define NIPOINTS 20

SplinePoint *TraceAndFitSpline(StrokeContext *c, Spline *s, bigreal t_start, bigreal t_end,
    SplinePoint *start, int nci_hint, int is_right) {
    NibOffset no;
    TPoint tp[NIPOINTS];
    bigreal nidiff, t, nt;
    int i, no_idx;
    BasePoint xy, ut, txy, ut_s, ut_e;
    SplinePoint *sp;

    nidiff = (t_end - t_start) / (NIPOINTS-1);

    for ( i=0, t=t_start; i<NIPOINTS; ++i, t+= nidiff ) {
	xy = SPLINEPVAL(s, t);
	ut = SplineUTanVecAt(s, t);
	if ( i==0 )
	    ut_s = ut;
	if ( i==(NIPOINTS-1) )
	    ut_e = ut;
	CalcNibOffset(c, BP_REV_IF(is_right, ut), &no, nci_hint);
	// We shouldn't hit a nib corner in the middle of tracing
	assert(i==0||i==(NIPOINTS-1)||no.nci[0]==no.nci[1]);
	// XXX change
	no_idx = 0;
	nci_hint = no.nci[no_idx];
	tp[i].x = xy.x + no.off[no_idx].x;
	tp[i].y = xy.y + no.off[no_idx].y;
	tp[i].t = (bigreal)i/(NIPOINTS-1);
	printf("TPoints i=%d, x=%lf, y=%lf, t=%lf\n", i, tp[i].x, tp[i].y, tp[i].t);
    }
    // XXX not exact enough assert(BPWITHIN(start->me, ((BasePoint) { tp[0].x, tp[0].y }), 1e-4));
    start->nextcp = BP_ADD(start->me, ut_s);
    start->nonextcp = false;
    sp = SplinePointCreate(tp[NIPOINTS-1].x, tp[NIPOINTS-1].y);
    sp->prevcp = BP_ADD(sp->me, BP_REV(ut_e));
    sp->noprevcp = false;
    ApproximateSplineFromPointsSlopes(start,sp,tp+1,NIPOINTS-2,false);

    return sp;
}

static bigreal NextNibT(StrokeContext *c, Spline *s, int nci, bigreal cur_t, BasePoint *cur_ut, int reverse) {
    int si[3], j, cnt = 2;
    bigreal t, min_t = 2;
    BasePoint ut = *cur_ut, min_ut;

    // Look for angles at higher and lower angle corners.
    si[0] = (nci+1)%c->n;
    si[1] = nci;
    // If we're on a corner slope we need to look for corners on either side
    // but also the current corner because of inflection points.
    if ( BPNEAR(c->nibcorners[nci].utanvec, BP_REV_IF(reverse, ut)) ) {
	si[2] = (c->n+nci-1)%c->n;
	++cnt;
    }
    printf("i: %d, si0: %d, si1: %d, si2: %d\n", nci, si[0], si[1], si[2]);

    for ( j=0; j<cnt; ++j ) {
	ut = c->nibcorners[si[j]].utanvec;
	ut = BP_REV_IF(reverse, ut);
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

static SplinePoint *AddJoint(StrokeContext *c, BasePoint sxy, SplinePoint *start_p,
                             NibOffset *nos, int nos_idx, NibOffset *noe, int noe_idx, int bk) {
    SplinePoint *sp;
    // assert( BPNEAR((BP_ADD(nos->off[nos_idx], sxy)), start_p->me) );

    printf("nos_idx=%d, noe_idx=%d, snci=%d, enci=%d\n", nos_idx, noe_idx, nos->nci[nos_idx], noe->nci[noe_idx]);
    sp = AppendCubicSplineSetPortion(c->nibcorners[nos->nci[nos_idx]].on_nib, nos->nt[nos_idx],
                                     c->nibcorners[noe->nci[noe_idx]].on_nib, noe->nt[noe_idx], start_p, bk);
    return sp;
}

/*
static void AddCap(StrokeContext *c, BasePoint sxy, SplineSet *ss, NibOffset *nos, NibOffset *noe) {
    SplinePoint *sp;
    NibOffset no1, no2;
    BasePoint sxy;
    int linear;

	    sxy = SPLINEPVAL(s, 0.0);
	    CalcNibOffset(c, BP_REV_IF(is_right, ut_start), &no, -1);
	    no_idx = cw_start ? NIBOFF_HIGH_IDX : NIBOFF_LOW_IDX;
	    oxy = BP_ADD(sxy, no.off[no_idx]);
	linear = SplineIsLinearish(s);
    CalcNibOffset(c, ut_end, &no1, -1);
    CalcNibOffset(c, BP_REV(ut_end), &no2, -1);
    return sp;
} */

static SplineSet *SplinesToContours(SplineSet *ss, StrokeContext *c) {
    NibOffset no, no_last;
    Spline *s, *first=NULL;
    SplineSet *left=NULL, *right=NULL, *cur;
    SplinePoint *sp;
    BasePoint ut_ini = { 2,0 }, ut_start, ut_mid, ut_end, ut_endlast, ut_diff;
    BasePoint oxy, sxy;
    bigreal last_t, t;
    int is_right, linear;
    int no_idx, nol_idx, closed = ss->first->prev!=NULL;

    if ( !c->remove_outer || !closed )
	left = chunkalloc(sizeof(SplineSet));
    if ( !c->remove_inner || !closed )
	right = chunkalloc(sizeof(SplineSet));

    for ( s=ss->first->next; s!=NULL && s!=first; s=s->to->next ) {
	if ( first==NULL )
	    first = s;

	if ( SplineLength(s)==0 )
	    continue;		/* We can safely ignore it because it is of zero length */

	ut_start = SplineUTanVecAt(s, 0.0);
	if ( ut_ini.x==2 )
	    ut_ini = ut_start;
	linear = SplineIsLinearish(s);
	if ( linear ) {
	    ut_end = ut_start;
	    no_idx = NIBOFF_HIGH_IDX;
	} else {
	    ut_end = SplineUTanVecAt(s, 1.0);
	    no_idx = SplineTurningCWAt(s, 0.0) ? NIBOFF_HIGH_IDX : NIBOFF_LOW_IDX;
	}

	// Left then right
	for ( is_right=0; is_right<=1; ++is_right ) {

	    if ( is_right ) {
		if ( right==NULL )
		    continue;
		cur = right;
	    } else {
		if ( left==NULL )
		    continue;
		cur = left;
	    }

	    sxy = SPLINEPVAL(s, 0.0);
	    CalcNibOffset(c, BP_REV_IF(is_right, ut_start), &no, -1);
	    oxy = BP_ADD(sxy, no.off[no_idx]);

	    // Create or verify initial spline point
	    printf("is_right=%d nci=%d nc.x=%lf nc.y=%lf x=%lf, y=%lf, ut_start=%lf,%lf\n", is_right, no.nci[no_idx], c->nibcorners[no.nci[no_idx]].on_nib->me.x, c->nibcorners[no.nci[no_idx]].on_nib->me.y, oxy.x, oxy.y, ut_start.x, ut_start.y);

	    if ( cur->first==NULL ) {
		cur->first = SplinePointCreate(oxy.x, oxy.y);
		cur->last = cur->first;
	    } else if ( !BPNEAR(cur->last->me, oxy) ) {
		ut_diff = UTanVecDiff(ut_endlast, ut_start);
		if ( (ut_diff.y > 0) == !is_right ) {
		    sp = SplinePointCreate(oxy.x, oxy.y);
		    SplineMake3(cur->last, sp);
		} else {
		    CalcNibOffset(c, ut_endlast, &no_last, -1);
		    sp = AddJoint(c, sxy, cur->last, &no_last, nol_idx, &no, no_idx, !is_right);
		}
		cur->last = sp;
	    }

	    // The path for this spline
	    if ( linear ) {
		sxy = SPLINEPVAL(s, 1.0);
		oxy = BP_ADD(sxy, no.off[no_idx]);
		sp = SplinePointCreate(oxy.x, oxy.y);
		SplineMake3(cur->last, sp);
		cur->last = sp;
	    } else {
		t = 0.0;
		ut_mid = ut_start;
		while ( t < 1.0 ) {
		    last_t = t;
		    t = NextNibT(c, s, no.nci[no_idx], t, &ut_mid, is_right);
		    assert( t==-1.0 || t > last_t );
		    if ( t==-1.0 ) // No relevant slope found -- copy rest of spline
			t = 1.0;
		    printf("nci: %d, ut_mid=%lf,%lf, t=%.15lf\n", no.nci[no_idx], ut_mid.x, ut_mid.y, t);
		    if ( no.linear )
			sp = AppendCubicSplinePortion(s, last_t, t, cur->last);
		    else
			sp = TraceAndFitSpline(c, s, last_t, t, cur->last, no.nci[no_idx], is_right);
		    cur->last = sp;
		    if (t < 1.0) {
			sxy = SPLINEPVAL(s, t);
			no_idx = SplineTurningCWAt(s, t) ? NIBOFF_HIGH_IDX : NIBOFF_LOW_IDX;
			CalcNibOffset(c, BP_REV_IF(is_right, ut_mid), &no, no.nci[NIBOFF_LOW_IDX]);
			printf("lowidx=%d, lowx=%lf, lowy=%lf, highidx=%d, highx=%lf, highy=%lf\n", no.nci[1], no.off[1].x, no.off[1].y, no.nci[0], no.off[0].x, no.off[0].y);
			oxy = BP_ADD(sxy, no.off[no_idx]);
			if ( !BPNEAR(cur->last->me, oxy) ) {
			    // One of the two points should be the same
			    //assert(BPNEAR(cur->last->me, BP_ADD(sxy, no.off[(no_idx+1)%2])));
			    sp = SplinePointCreate(oxy.x, oxy.y);
			    SplineMake3(cur->last, sp);
			    cur->last = sp;
			}
		    }
		}
	    }
	}
	ut_endlast = ut_end;
        nol_idx = SplineTurningCWAt(s, 1.0) ? NIBOFF_HIGH_IDX : NIBOFF_LOW_IDX;
    }
    if ( (left!=NULL && left->first==NULL) || (right!=NULL && right->first==NULL) ) {
	// Presumably the path had only zero-length splines
	assert( (left==NULL || left->first==NULL) && (right==NULL || right->first==NULL) );
	chunkfree(left, sizeof(SplineSet));
	chunkfree(right, sizeof(SplineSet));
	return NULL;
    }
    if ( !closed ) {
	SplineSetReverse(right);
/*	AddCap(c, ss, left, false);
	left->next = right;
	right = NULL;
	SplineSetJoin(left, true, 1e-8, &closed);
	assert(closed);
	AddCap(c, ss, left, true);
        CalcNibOffset(c, ut_start, &no, -1);
        CalcNibOffset(c, BP_REV(ut_start), &no2, -1);
	AddCap(c, left->last, &no, &no2);
	SplineSetJoin(left, true, 1e-8, &closed);
	assert(closed); */
	SplineMake3(left->last, right->first);
	SplineMake3(right->last, left->first);
	right = NULL;
	left->last = left->first;
    } else {
	if ( left!=NULL ) {
            left = SplineSetJoin(left, true, 1e-8, &closed);
	    assert(closed);
	}
	if ( right!=NULL ) {
            right = SplineSetJoin(right, true, 1e-8, &closed);
	    assert(closed);
	    if ( left != NULL ) {
		SplineSetReverse(right);
		left->next = right;
	    } else
		left = right;
	    right = NULL;
	}
    }
    return left;
}

/******************************************************************************/
/* ******************************* Unit Stuff ******************************* */
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

/******************************************************************************/
/* ******************************* Higher-Level ***************************** */
/******************************************************************************/

static SplinePointList *SinglePointStroke(SplinePoint *sp,struct strokecontext *c) {
    SplineSet *ret;
    SplinePoint *sp1, *sp2;
    int i;
    real trans[6];

    if ( c->pentype==pt_circle && c->cap==lc_butt ) {
	ret = chunkalloc(sizeof(SplineSet));
	/* Leave as a single point */
	ret->first = ret->last = SplinePointCreate(sp->me.x,sp->me.y);
	ret->first->pointtype = pt_corner;
    } else { 
	memset(&trans, 0, sizeof(trans));
	trans[0] = trans[3] = 1;
	trans[4] = sp->me.x;
	trans[5] = sp->me.y;
	ret = SplinePointListCopy(c->nib);
	SplinePointListTransform(ret, trans, tpt_AllPoints);
    }
    return( ret );
}

static SplineSet *SplineSet_Stroke(SplineSet *ss,struct strokecontext *c, int order2) {
    SplineSet *base = ss;
    SplineSet *ret;
    int copied = 0;

    if ( base->first->next==NULL )
	ret = SinglePointStroke(base->first, c);
    else {
	if ( base->first->next->order2 ) {
	    base = SSPSApprox(ss);
	    copied = 1;
	}
	ret = SplinesToContours(base, c);
    }
    if ( order2 )
	ret = SplineSetsConvertOrder(ret,order2);
    if ( copied )
	SplinePointListFree(base);
    return(ret);
}

static SplineSet *SplineSets_Stroke(SplineSet *ss,struct strokecontext *c, int order2) {
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

static int BuildNibCorners(StrokeContext *c, int maxp) {
    int i, max_utan_index = -1;
    BasePoint max_utanangle = { -1, -DBL_MIN };
    SplinePoint *sp;
    NibCorner *tpc;

    if ( c->nibcorners==NULL )
	c->nibcorners = calloc(maxp,sizeof(NibCorner));

    for ( sp=c->nib->first, i=0; ; ) {
	if ( i==maxp ) { // We guessed wrong
	    maxp *= 2;
	    tpc = calloc(maxp, sizeof(NibCorner));
	    memcpy(tpc, c->nibcorners, (i-1)*sizeof(NibCorner));
	    free(c->nibcorners);
	    c->nibcorners = tpc;
	}
	c->nibcorners[i].on_nib = sp;
	c->nibcorners[i].linear = SplineIsLinear(sp->prev);
	c->nibcorners[i].utanvec = SplineUTanVecAt(sp->prev, 0.0);
	if ( UTanAngleGreater(c->nibcorners[i].utanvec, max_utanangle) ) {
	    max_utan_index = i;
	    max_utanangle = c->nibcorners[i].utanvec;
	}
	++i;
	sp=sp->next->to;
	if ( sp==c->nib->first )
	    break;
    }
    c->n = i;

    if (max_utan_index != 0) {
	tpc = malloc(c->n*sizeof(NibCorner));

	memcpy(tpc, c->nibcorners, c->n*sizeof(NibCorner));
	memcpy(tpc, c->nibcorners+max_utan_index, (c->n-max_utan_index)*sizeof(NibCorner));
	memcpy(tpc+(c->n-max_utan_index), c->nibcorners, max_utan_index*sizeof(NibCorner));
	free(c->nibcorners);
	c->nibcorners = tpc;
    }

    for (i = 0; i < c->n; i++)
        printf("Corner %d: x=%lf, y=%lf, ut=%lf,%lf\n", i, c->nibcorners[i].on_nib->me.x,
	       c->nibcorners[i].on_nib->me.y, c->nibcorners[i].utanvec.x, c->nibcorners[i].utanvec.y);

    return maxp;
}

SplineSet *SplineSetStroke(SplineSet *ss,StrokeInfo *si, int order2) {
    StrokeContext c;
    SplineSet *nibs, *nib, *bnext, *first, *last, *cur;
    bigreal sn, co;
    DBounds b;
    real trans[6];
    int max_pc;

    if ( si->stroke_type==si_centerline )
	IError("centerline not handled");

    memset(&c,0,sizeof(c));
    c.join = si->join;
    c.cap  = si->cap;
    c.miterlimit = /* -cos(theta) */ -.98 /* theta=~11 degrees, PS miterlimit=10 */;
    c.radius = si->radius;
    c.remove_inner = si->removeinternal;
    c.remove_outer = si->removeexternal;
    c.leave_users_center = si->leave_users_center;
    c.scaled_or_rotated = si->factor!=NULL;

    c.leave_users_center = 1;

    memset(trans,0,sizeof(trans));
    if ( si->stroke_type==si_std || si->stroke_type==si_caligraphic ) {
	c.pentype = si->stroke_type==si_std ? pt_circle : pt_square;
	max_pc = 4;
	printf("Radius: %lf, Minor Radius: %lf\n", si->radius, si->minorradius);
	sn = sin(si->penangle);
	co = cos(si->penangle);
	trans[0] = si->radius*co;
	trans[1] = si->radius*sn;
	trans[2] = si->minorradius*-sn;
	trans[3] = si->minorradius*co;
	nibs = UnitShape(si->stroke_type==si_std ? 0 : -4);
	SplinePointListTransformExtended(nibs,trans,tpt_AllPoints,tpmask_dontTrimValues);
    } else {
	c.pentype = pt_poly;
	max_pc = 20; // a guess
	nibs = si->nib;
	if ( !c.leave_users_center ) {
	    SplineSetQuickBounds(nibs,&b);
	    trans[0] = trans[3] = 1;
	    trans[4] = -(b.minx+b.maxx)/2;
	    trans[5] = -(b.miny+b.maxy)/2;
	    SplinePointListTransformExtended(nibs,trans,tpt_AllPoints,tpmask_dontTrimValues);
	}
    }

    memset(trans,0,sizeof(trans));
    trans[0] = trans[3] = 1;
    first = last = NULL;
    for ( nib=nibs; nib!=NULL; nib=nib->next ) {
	int reversed = false;
	if ( !SplinePointListIsClockwise(nib) ) {
	    reversed = true;
	    SplineSetReverse(nib);
	}
	c.nib = nib;
	max_pc = BuildNibCorners(&c, max_pc);
	cur = SplineSets_Stroke(ss,&c,order2);
	if ( reversed ) {
	    SplineSet *ss;
	    for ( ss=cur; ss!=NULL; ss=ss->next )
		SplineSetReverse(ss);
	    SplineSetReverse(nib);
	}
        if ( first==NULL )
	    first = last = cur;
	else
	    last->next = cur;
	if ( last!=NULL )
	    while ( last->next!=NULL )
	        last = last->next;
    }

    free(c.nibcorners);

    if ( si->stroke_type==si_std || si->stroke_type==si_caligraphic ) {
	SplinePointListFree(si->nib);
	si->nib = NULL;
    }

    return( first );
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
