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

// Rounding makes even slope-continuous splines only approximately so. 
#define INTERSPLINE_MARGIN .5
#define INTRASPLINE_MARGIN (1e-8)
#define FIXUP_MARGIN (1e-2)

#define BP_LENGTHSQ(v) ((v).x*(v).x+(v).y*(v).y)
#define BP_REV(v) (BasePoint) { -(v).x, -(v).y }
#define BP_REV_IF(t, v) (t ? (BasePoint) { -(v).x, -(v).y } : (v))
#define BP_ADD(v1, v2) (BasePoint) { (v1).x + (v2).x, (v1).y + (v2).y } 
#define BP_UNINIT ((BasePoint) { -INFINITY, INFINITY })
#define BP_IS_UNINIT(v) ((v).x==-INFINITY && (v).y==INFINITY)
#define NORMANGLE(a) ((a)>PI?(a)-2*PI:(a)<-PI?(a)+2*PI:(a))

#define UTMIN ((BasePoint) { -1, -DBL_MIN })
#define BPNEAR(bp1, bp2) BPWITHIN(bp1, bp2, INTRASPLINE_MARGIN)

enum pentype { pt_circle, pt_square, pt_convex };

/* c->nibcorners is a structure that models the "corners" of the current nib
 * treated as if all splines were linear, and therefore as a convex polygon.
 */

#define NC_IN_IDX 0
#define NC_OUT_IDX 1

#define NC_NEXTI(c, i) ((i+1)%c->n)
#define NC_PREVI(c, i) ((c->n+i-1)%c->n)

typedef struct nibcorner {
    SplinePoint *on_nib;
    BasePoint utv[2];    // Unit tangents entering and leaving corner
    unsigned int linear: 1;
} NibCorner;

typedef struct strokecontext {
    enum pentype pentype;
    enum linejoin join;
    enum linecap cap;
    enum stroke_rmov rmov;
    bigreal miterlimit;	// For miter joins
	    // PostScript uses 1/sin( theta/2 ) as their miterlimit
	    //  I use -cos(theta). (where theta is the angle between the slopes
	    //  same idea, different implementation
    SplineSet *nib;
    int n;
    NibCorner *nibcorners;
    unsigned int remove_inner: 1;
    unsigned int remove_outer: 1;
    unsigned int leave_users_center: 1;	// Don't center the nib when stroking
    unsigned int simplify:1;
    unsigned int extrema:1;
} StrokeContext;

#define NIBOFF_CW_IDX 0
#define NIBOFF_CCW_IDX 1

typedef struct niboffset {
    BasePoint utanvec;        // The reference angle
    int nci[2];
    BasePoint off[2];
    bigreal nt;               // The t value on the nib spline for the angle
    unsigned int at_line: 1;  // Whether the angle is ambiguous
    unsigned int curve: 1;    // Whether the angle is on a curve or line
    unsigned int reversed: 1; // Whether the calculations were reversed for
                              // tracing the right side
} NibOffset;

enum jointype { jt_flat, jt_joint, jt_cap };

static char *glyphname=NULL;

/* Orders unit tangent vectors on an absolute -PI+e to PI
 * equivalent basis
 */
static int UTanVecGreater(BasePoint uta, BasePoint utb) {
    assert(    RealNear(BP_LENGTHSQ(uta),1)
            && RealNear(BP_LENGTHSQ(utb),1) );

    if (uta.y >= 0) {
	if (utb.y < 0)
	    return true;
	return uta.x < utb.x && !BPNEAR(uta, utb);
    }
    if (utb.y >= 0)
	return false;
    return uta.x > utb.x && !BPNEAR(uta, utb);
}

/* True if rotating from ut1 to ut3 in the specified direction
 * the angle passes through ut2.
 *
 * Note: Returns true if ut1 is near ut2 but false if ut2 is
 * near ut3
 */
static int UTanVecsSequent(BasePoint ut1, BasePoint ut2, BasePoint ut3,
                           int ccw) {
    BasePoint tmp;

    if ( BPNEAR(ut1, ut2) )
	return true;

    if ( BPNEAR(ut2, ut3) || BPNEAR(ut1, ut3) )
	return false;

    if (ccw) {
	tmp = ut1;
	ut1 = ut3;
	ut3 = tmp;
    }

    if ( UTanVecGreater(ut1, ut3) ) {
	return UTanVecGreater(ut1, ut2) && UTanVecGreater(ut2, ut3);
    } else {
	return    (UTanVecGreater(ut1, ut2) && UTanVecGreater(ut3, ut2))
	       || (UTanVecGreater(ut2, ut1) && UTanVecGreater(ut2, ut3));
    }
}

static int JointBendsCW(BasePoint ut_ref, BasePoint ut_vec) {
    // magnitude of cross product
    bigreal r = ut_ref.x*ut_vec.y - ut_ref.y*ut_vec.x;
    if ( RealWithin(r, 0.0, 1e-8) )
	return false;
    return r > 0;
}

static int BuildNibCorners(StrokeContext *c, int maxp) {
    int i, max_utan_index = -1;
    BasePoint max_utanangle = UTMIN;
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
	c->nibcorners[i].linear = SplineIsLinear(sp->next);
	c->nibcorners[i].utv[NC_IN_IDX] = SplineUTanVecAt(sp->prev, 1.0);
	c->nibcorners[i].utv[NC_OUT_IDX] = SplineUTanVecAt(sp->next, 0.0);
	if ( UTanVecGreater(c->nibcorners[i].utv[NC_IN_IDX], max_utanangle) ) {
	    max_utan_index = i;
	    max_utanangle = c->nibcorners[i].utv[NC_IN_IDX];
	}
	++i;
	sp=sp->next->to;
	if ( sp==c->nib->first )
	    break;
    }
    c->n = i;

    // Put in order of decreasing utanvec[NC_IN_IDX]
    if (max_utan_index != 0) {
	tpc = malloc(c->n*sizeof(NibCorner));

	memcpy(tpc, c->nibcorners, c->n*sizeof(NibCorner));
	memcpy(tpc, c->nibcorners+max_utan_index,
	       (c->n-max_utan_index)*sizeof(NibCorner));
	memcpy(tpc+(c->n-max_utan_index), c->nibcorners,
	       max_utan_index*sizeof(NibCorner));
	free(c->nibcorners);
	c->nibcorners = tpc;
    }

    for (i = 0; i < c->n; i++)
        printf("Corner %d: x=%lf, y=%lf, utin=%lf,%lf, utout=%lf, %lf\n",
	       i, c->nibcorners[i].on_nib->me.x,
	       c->nibcorners[i].on_nib->me.y,
	       c->nibcorners[i].utv[NC_IN_IDX].x,
	       c->nibcorners[i].utv[NC_IN_IDX].y,
	       c->nibcorners[i].utv[NC_OUT_IDX].x,
	       c->nibcorners[i].utv[NC_OUT_IDX].y);

    return maxp;
}

// The index into c->nibcorners associated with unit tangent ut
static int IndexForUTanVec(StrokeContext *c, BasePoint ut, int nci_hint) {
    int nci;

    assert( nci_hint == -1 || (nci_hint >= 0 && nci_hint < c->n ) );

    if (   nci_hint != -1
        && UTanVecsSequent(c->nibcorners[nci_hint].utv[NC_IN_IDX], ut,
	                   c->nibcorners[NC_NEXTI(c, nci_hint)].utv[NC_IN_IDX],
                           false) )
	nci = nci_hint;
    else // XXX Replace with binary search
	for (nci = 0; nci<c->n; ++nci)
	    if ( UTanVecsSequent(c->nibcorners[nci].utv[NC_IN_IDX], ut,
	                         c->nibcorners[NC_NEXTI(c, nci)].utv[NC_IN_IDX],
	                         false) )
		break;

    assert( UTanVecsSequent(c->nibcorners[nci].utv[NC_IN_IDX], ut,
                            c->nibcorners[NC_NEXTI(c, nci)].utv[NC_IN_IDX],
                            false) );

    return nci;
}

/* The offset from the stroked curve associated with unit tangent ut
 * (or it's reverse).
 */
static NibOffset *CalcNibOffset(StrokeContext *c, BasePoint ut, int reverse,
                                NibOffset *no, int nci_hint) {
    int nci, ncni, ncpi;
    Spline *ns;
    BasePoint tmp;

    if ( no==NULL )
	no = malloc(sizeof(NibOffset));

    memset(no,0,sizeof(NibOffset));
    no->utanvec = ut; // Store the unreversed value for reference

    if ( reverse ) {
	ut = BP_REV(ut);
	no->reversed = 1;
    }

    nci = no->nci[0] = no->nci[1] = IndexForUTanVec(c, ut, nci_hint);
    ncpi = NC_PREVI(c, nci);
    ncni = NC_NEXTI(c, nci);

    // The only case where one of two points might draw the same angle,
    // and therefore the only case where the array values differ
    if (   c->nibcorners[ncpi].linear
        && BPNEAR(ut, c->nibcorners[ncpi].utv[NC_OUT_IDX]) ) {
	no->nt = 0.0;
	no->off[NIBOFF_CCW_IDX] = c->nibcorners[ncpi].on_nib->me;
	no->nci[NIBOFF_CCW_IDX] = ncpi;
	no->off[NIBOFF_CW_IDX] = c->nibcorners[nci].on_nib->me;
	no->at_line = true;
	no->curve = false;
    // Whether ut is (effectively) between IN and OUT, in 
    // which case the point draws the offset curve
    } else if (   UTanVecsSequent(c->nibcorners[nci].utv[NC_IN_IDX], ut,
                                  c->nibcorners[nci].utv[NC_OUT_IDX], false) 
               || BPNEAR(ut, c->nibcorners[nci].utv[NC_OUT_IDX] ) ) {
	no->nt = 0.0;
	no->off[0] = no->off[1] = c->nibcorners[nci].on_nib->me;
	no->curve = false;
    } else {
	// ut is on a spline curve clockwise from point nci
	assert( UTanVecsSequent(c->nibcorners[nci].utv[NC_OUT_IDX], ut,
	                        c->nibcorners[ncni].utv[NC_IN_IDX], false) );
	// Nib splines are locally convex and therefore have t value per slope
	ns = c->nibcorners[nci].on_nib->next;
	no->nt = SplineSolveForUTanVec(ns, ut, 0.0);
	no->off[0] = no->off[1] = SPLINEPVAL(ns, no->nt);
	no->curve = true;
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
    Spline *s;
    int n = 1;
    bigreal anglesum = 0, angle, pangle, d, pcpangle, ncpangle;

    if ( ss->first==NULL )
	return Shape_TooFewPoints;
    if ( ss->first->prev==NULL )
	return Shape_NotClosed;
    if ( !SplinePointListIsClockwise(ss) )
	return Shape_CCW;

    s = ss->first->prev;
    pangle = atan2(s->to->me.y - s->from->me.y, s->to->me.x - s->from->me.x);
    s = ss->first->next;
    SplinePointListSelect(ss, false);
    while ( true ) {
	s->from->selected = true;
	if ( BPWITHIN(s->from->me, s->to->me,1e-2) )
	    return Shape_TinySpline;
	angle = atan2(s->to->me.y - s->from->me.y, s->to->me.x - s->from->me.x);
	if ( RealWithin(angle, pangle, 1e-4) )
	    return Shape_PointOnEdge;
	d = pangle-angle;
	d = NORMANGLE(d);
	if ( d<0 )
	    return Shape_CCWTurn;
	anglesum += d;
	if ( !s->from->noprevcp ) {
	    pcpangle = atan2(s->to->prevcp.y - s->to->me.y,
	                     s->to->prevcp.x - s->to->me.x);
	    d = pcpangle-pangle;
	    d = NORMANGLE(d);
	    printf("prevcp angle: %lf, d1: %lf, d2:%lf\n", pcpangle, d, NORMANGLE(pcpangle-angle));
	} else
	    pcpangle = pangle;
	if ( !s->from->nonextcp ) {
	    ncpangle = atan2(s->to->nextcp.y - s->to->me.y,
	                     s->to->nextcp.x - s->to->me.x);
	    d = ncpangle-pcpangle;
	    d = NORMANGLE(d);
	    printf("nextcp angle: %lf, d1: %lf, d2:%lf\n", ncpangle, d, NORMANGLE(angle-ncpangle));
	}
	s->from->selected = false;
	s=s->to->next;
	pangle = angle;
	if ( s==ss->first->next )
	    break;
	++n;
    }
    if ( n<3 )
	return Shape_TooFewPoints;
    if ( !RealWithin(anglesum, 2*PI, 1e-1) )
	return Shape_SelfIntersects;
    return Shape_Convex;
}

/* Copies the portion of s from t_start to t_end and then translates
 * and appends it to t_start. The new end point is returned. "Reversing"
 * t_start and t_end reverses the copy's direction.
 *
 * Calculations cribbed from https://stackoverflow.com/a/879213
 */
SplinePoint *AppendCubicSplinePortion(Spline *s, bigreal t_start, bigreal t_end,
                                 SplinePoint *start) {
    extended u_start = 1-t_start, u_end = 1-t_end;
    SplinePoint *end;
    BasePoint v, qf, qcf, qct, qt;

    // XXX Maybe this should be length based
    if ( RealWithin(t_start, t_end, 1e-3) )
	return start;

    // Intermediate calculations
    qf.x =    s->from->me.x*u_start*u_start
            + s->from->nextcp.x*2*t_start*u_start
            + s->to->prevcp.x*t_start*t_start;
    qcf.x =   s->from->me.x*u_end*u_end
            + s->from->nextcp.x*2*t_end*u_end
            + s->to->prevcp.x*t_end*t_end;
    qct.x =   s->from->nextcp.x*u_start*u_start
            + s->to->prevcp.x*2*t_start*u_start
            + s->to->me.x*t_start*t_start;
    qt.x =    s->from->nextcp.x*u_end*u_end
            + s->to->prevcp.x*2*t_end*u_end
            + s->to->me.x*t_end*t_end;

    qf.y =    s->from->me.y*u_start*u_start
	    + s->from->nextcp.y*2*t_start*u_start
            + s->to->prevcp.y*t_start*t_start;
    qcf.y =   s->from->me.y*u_end*u_end
            + s->from->nextcp.y*2*t_end*u_end
            + s->to->prevcp.y*t_end*t_end;
    qct.y =   s->from->nextcp.y*u_start*u_start
            + s->to->prevcp.y*2*t_start*u_start
            + s->to->me.y*t_start*t_start;
    qt.y =    s->from->nextcp.y*u_end*u_end
            + s->to->prevcp.y*2*t_end*u_end
            + s->to->me.y*t_end*t_end;

    // Difference vector to offset other points
    v.x = start->me.x - (qf.x*u_start + qct.x*t_start);
    v.y = start->me.y - (qf.y*u_start + qct.y*t_start);
    //printf("vx = %lf, vy = %lf\n", v.x, v.y);

    end = SplinePointCreate(qcf.x*u_end + qt.x*t_end + v.x,
                            qcf.y*u_end + qt.y*t_end + v.y);

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

/* Copies the splines between s_start (at t_start) and s_end (at t_end)
 * after point dst_start. Backward (and therefore forward) is relative
 * to next/prev rather than clockwise/counter-clockwise.
 */
SplinePoint *AppendCubicSplineSetPortion(Spline *s_start, bigreal t_start,
                                         Spline *s_end, bigreal t_end,
					 SplinePoint *dst_start, int backward) {
    Spline *s;

    if ( backward && RealWithin(t_start, 0.0, 1e-4) ) {
	t_start = 1;
	s_start = s_start->from->prev;
    } else if ( !backward && RealWithin(t_start, 1.0, 1e-4) ) {
	t_start = 0;
	s_start = s_start->to->next;
    }
    if ( backward && RealWithin(t_end, 1.0, 1e-4) ) {
	t_end = 0.0;
	s_end = s_end->to->next;
    } else if ( !backward && RealWithin(t_end, 0.0, 1e-4) ) {
	t_end = 1.0;
	s_end = s_end->from->prev;
    }
    s = s_start;

    // Handle the single spline case
    if (    s==s_end
         && (( t_start<t_end && !backward ) || (t_start>t_end && backward))) {
	dst_start = AppendCubicSplinePortion(s, t_start, t_end, dst_start);
	return dst_start;
    }

    dst_start = AppendCubicSplinePortion(s, t_start, backward ? 0 : 1,
                                         dst_start);

    while ( 1 ) {
        s = backward ? s->from->prev : s->to->next;
	if ( s==s_end )
	    break;
	assert( s!=NULL && s!=s_start ); // XXX turn into runtime check
        dst_start = AppendCubicSplinePortion(s, backward ? 1 : 0,
	                                     backward ? 0 : 1, dst_start);
    }
    dst_start = AppendCubicSplinePortion(s, backward ? 1 : 0, t_end, dst_start);
    return dst_start;
}

#define NIPOINTS 30

SplinePoint *TraceAndFitSpline(StrokeContext *c, Spline *s, bigreal t_start,
                               bigreal t_end, SplinePoint *start,
			       int nci_hint, int is_right) {
    NibOffset no;
    TPoint tp[NIPOINTS];
    bigreal nidiff, t, nt;
    int i, is_ccw;
    BasePoint xy, ut, txy, ut_s, ut_e;
    SplinePoint *sp;

    nidiff = (t_end - t_start) / (NIPOINTS-1);

    is_ccw = SplineTurningCCWAt(s, t_start);
    no.nci[0] = no.nci[1] = nci_hint;
    for ( i=0, t=t_start; i<NIPOINTS; ++i, t+= nidiff ) {
	xy = SPLINEPVAL(s, t);
	ut = SplineUTanVecAt(s, t);
	if ( i==0 )
	    ut_s = ut;
	if ( i==(NIPOINTS-1) )
	    ut_e = ut;
	CalcNibOffset(c, ut, is_right, &no, no.nci[is_ccw]);
	tp[i].x = xy.x + no.off[is_ccw].x;
	tp[i].y = xy.y + no.off[is_ccw].y;
	tp[i].t = (bigreal)i/(NIPOINTS-1);
	// printf("TPoints i=%d, t(x=%lf, y=%lf, tp.t=%lf), spline t=%lf sxy=%lf,%lf, off=%lf,%lf\n, ut=%lf,%lf", i, tp[i].x, tp[i].y, tp[i].t, t, xy.x, xy.y, no.off[is_ccw].x, no.off[is_ccw].y, ut.x, ut.y);
    }
    if ( !BPWITHIN(start->me, ((BasePoint) { tp[0].x, tp[0].y }),
                   FIXUP_MARGIN) )
	LogError(_("Warning: Coordinate diff %lf greater than margin %lf\n"),
                 fmax(fabs(start->me.x-tp[0].x), fabs(start->me.y-tp[0].y)),
		 FIXUP_MARGIN);
    start->nextcp = BP_ADD(start->me, ut_s);
    start->nonextcp = false;
    sp = SplinePointCreate(tp[NIPOINTS-1].x, tp[NIPOINTS-1].y);
    sp->prevcp = BP_ADD(sp->me, BP_REV(ut_e));
    sp->noprevcp = false;
    ApproximateSplineFromPointsSlopes(start,sp,tp+1,NIPOINTS-2,false);

    return sp;
}

static BasePoint SplineStrokeNextAngle(StrokeContext *c, BasePoint ut,
                                       int is_ccw, int *curved, int reverse,
				       int nci_hint) {
    int nci, ncni, ncpi, inout;

    ut = BP_REV_IF(reverse, ut);
    nci = IndexForUTanVec(c, ut, nci_hint);
    ncni = NC_NEXTI(c, nci);
    ncpi = NC_PREVI(c, nci);

    if ( BPNEAR(ut, c->nibcorners[nci].utv[NC_IN_IDX]) ) {
	if ( is_ccw ) {
	    if ( BPNEAR(c->nibcorners[nci].utv[NC_IN_IDX],
	                c->nibcorners[ncpi].utv[NC_OUT_IDX]) ) {
		assert(c->nibcorners[ncpi].linear);
		inout = NC_IN_IDX;
		*curved = false;
	    } else {
		inout = NC_OUT_IDX;
		*curved = true;
	    }
	    nci = ncpi;
	} else { // CW
	    if ( BPNEAR(c->nibcorners[nci].utv[NC_IN_IDX],
	                c->nibcorners[nci].utv[NC_OUT_IDX]) ) {
		nci = ncni;
		inout = NC_IN_IDX;
		*curved = true;
	    } else {
		inout = NC_OUT_IDX;
		*curved = false;
	    }
	}
    } else if ( BPNEAR(ut, c->nibcorners[nci].utv[NC_OUT_IDX]) ) {
	if ( is_ccw ) {
	    if ( BPNEAR(c->nibcorners[nci].utv[NC_IN_IDX],
	                c->nibcorners[nci].utv[NC_OUT_IDX]) ) {
		nci = ncpi;
		inout = NC_OUT_IDX;
		*curved = true;
	    } else {
		inout = NC_IN_IDX;
		*curved = false;
	    }
	} else { // CW
	    if ( BPNEAR(c->nibcorners[nci].utv[NC_OUT_IDX],
	                c->nibcorners[ncni].utv[NC_IN_IDX]) ) {
		assert(c->nibcorners[nci].linear);
		inout = NC_OUT_IDX;
		*curved = false;
	    } else {
		inout = NC_IN_IDX;
		*curved = true;
	    }
	    nci = ncni;
	}
    } else if ( UTanVecsSequent(c->nibcorners[nci].utv[NC_IN_IDX], ut,
                c->nibcorners[nci].utv[NC_OUT_IDX], false) ) {
	if ( is_ccw )
	    inout = NC_IN_IDX;
	else
	    inout = NC_OUT_IDX;
	*curved = false;
    } else {
	assert( UTanVecsSequent(c->nibcorners[nci].utv[NC_OUT_IDX], ut,
	        c->nibcorners[ncni].utv[NC_IN_IDX], false) );
	if ( is_ccw )
	    inout = NC_OUT_IDX;
	else {
	    nci = ncni;
	    inout = NC_IN_IDX;
	}
	*curved = true;
    }
    return BP_REV_IF(reverse, c->nibcorners[nci].utv[inout]);
}

static bigreal SplineStrokeNextT(StrokeContext *c, Spline *s, bigreal cur_t,
		                 int is_ccw, BasePoint *cur_ut,
				 int *curved, int reverse, int nci_hint) {
    int next_curved, nci, inout, at_line, icnt, i;
    bigreal next_t;
    extended poi[2], tp;
    BasePoint next_ut;

    assert( cur_ut!=NULL && curved!=NULL );

    next_ut = SplineStrokeNextAngle(c, *cur_ut, is_ccw, &next_curved,
                                    reverse, nci_hint);
    next_t = SplineSolveForUTanVec(s, next_ut, cur_t);

    // If there is an inflection point before next_t stop there first
    // so that tracing and splicing doesn't have to worry about inflections
    // and the traced curve will have fewer of them. The curved value returned
    // by SplineStrokeNextAngle is accurate for that portion of the curve
    //
    // The alternative would be to call
    //
    // SplineStrokeNextT(c, s, inflect_t, !is_ccw, cur_ut,
    //                   curved, reverse, nci_hint);
    //
    // to look for the next angle after the inflection turning the other way.
    if ( icnt = Spline2DFindPointsOfInflection(s, poi) ) {
	assert ( icnt < 2 || poi[0] < poi[1] );
	for ( i=0; i<2; ++i )
	    if (    poi[i] > cur_t
	         && !RealNear(poi[i], cur_t)
	         && (next_t==-1 || poi[i] < next_t) ) {
		next_t = poi[i];
		next_ut = SplineUTanVecAt(s, next_t);
		break;
	    }
    }

    printf("next_t:%lf, cur_t:%lf, next_curved:%d\n", next_t, cur_t,
           next_curved);

    if ( next_t==-1 ) {
	next_t = 1.0;
	*cur_ut = SplineUTanVecAt(s, next_t);
    } else
	*cur_ut = next_ut;

    *curved = next_curved;
    return next_t;
}

static SplinePoint *AddJoint(StrokeContext *c, SplinePoint *start_p,
                             NibOffset *nos, int is_ccw_s, NibOffset *noe,
			     int is_ccw_e, int bk) {
    SplinePoint *sp, *start_np, *end_np;

    start_np = c->nibcorners[nos->nci[is_ccw_s]].on_nib;
    end_np = c->nibcorners[noe->nci[is_ccw_e]].on_nib;
    sp = AppendCubicSplineSetPortion(start_np->next, nos->nt, end_np->next,
                                     noe->nt, start_p, bk);
    return sp;
}

/* Put the new endpoint exactly where the NibOffset calculation says it
 * should be to avoid cumulative append errors.
 */
static void SplineStrokeAppendFixup(SplinePoint *endp, BasePoint sxy,
                                    NibOffset *noe) {
    int i = 0;
    bigreal mg;
    BasePoint oxy[2], dxy[2];

    oxy[0] = BP_ADD(sxy, noe->off[0]);
    dxy[0] = BP_ADD(oxy[0], BP_REV(endp->me));

    // Use the closer point
    if ( noe->nci[0]!=noe->nci[1] ) {
	oxy[1] = BP_ADD(sxy, noe->off[1]);
	dxy[1] = BP_ADD(oxy[1], BP_REV(endp->me));
	if ( BP_LENGTHSQ(dxy[0]) > BP_LENGTHSQ(dxy[1]) )
	    i = 1;
    }

    mg = fmax(fabs(dxy[i].x), fabs(dxy[i].y));
    if ( mg > FIXUP_MARGIN ) {
	LogError(_("Warning: Coordinate diff %lf greater than margin %lf\n"),
	         mg, FIXUP_MARGIN);
	return;
    }

    endp->nextcp = BP_ADD(endp->nextcp, dxy[i]);
    endp->me = oxy[i];
}

static void HandleJoint(StrokeContext *c, SplineSet *cur,
                        BasePoint sxy, BasePoint oxy,
                        NibOffset *noi, int is_ccw,
                        BasePoint ut_endlast, int was_ccw,
                        int is_right) {
    NibOffset now;
    SplinePoint *sp;
    int jcw = JointBendsCW(ut_endlast, noi->utanvec);

    if ( BPWITHIN(cur->last->me, oxy, INTERSPLINE_MARGIN) ) {
	SplineStrokeAppendFixup(cur->last, sxy, noi);
    } else {
	if ( BPNEAR(ut_endlast, noi->utanvec) ) {
	    sp = SplinePointCreate(oxy.x, oxy.y);
	    SplineMake3(cur->last, sp);
	} else if ( jcw == !is_right ) {
	    // XXX Mark for cleanup
	    sp = SplinePointCreate(oxy.x, oxy.y);
	    SplineMake3(cur->last, sp);
	} else {
	    CalcNibOffset(c, ut_endlast, is_right, &now, -1);
	    sp = AddJoint(c, cur->last, &now, was_ccw, noi, is_ccw, jcw);
	    SplineStrokeAppendFixup(sp, sxy, noi);
	}
	cur->last = sp;
    }
}

static SplineSet *SplineContourOuterCCWRemoveOverlap(SplineSet *ss) {
    DBounds b;
    SplinePoint *sp;
    SplineSet *ss_tmp, *ss_last = NULL;

    SplineSetQuickBounds(ss,&b);
    ss_tmp = chunkalloc(sizeof(SplineSet));
    ss_tmp->first = SplinePointCreate(b.minx-100,b.miny-100);
    ss_tmp->last = SplinePointCreate(b.minx-100,b.maxy+100);
    SplineMake3(ss_tmp->first, ss_tmp->last);
    sp = SplinePointCreate(b.maxx+100,b.maxy+100);
    SplineMake3(ss_tmp->last, sp);
    ss_tmp->last = sp;
    sp = SplinePointCreate(b.maxx+100,b.miny-100);
    SplineMake3(ss_tmp->last, sp);
    SplineMake3(sp, ss_tmp->first);
    ss_tmp->last = ss_tmp->first;
    ss->next = ss_tmp;
    ss=SplineSetRemoveOverlap(NULL,ss,over_remove);
    for ( ss_tmp=ss; ss_tmp!=NULL; ss_last=ss_tmp, ss_tmp=ss_tmp->next )
	if (   ss_tmp->first->me.x==(b.minx-100)
	    || ss_tmp->first->me.x==(b.maxx+100) ) {
	    if ( ss_last==NULL )
		ss = ss->next;
    	    else
		ss_last->next = ss_tmp->next;
	    ss_tmp->next = NULL;
	    SplinePointListFree(ss_tmp);
	    return ss;
	}
    assert(0);
    return ss;
}

/* The guts of the stroking algorithm.
 */
static SplineSet *SplinesToContours(SplineSet *ss, StrokeContext *c) {
    NibOffset no, no_last;
    Spline *s, *first=NULL;
    SplineSet *left=NULL, *right=NULL, *cur;
    SplinePoint *sp;
    BasePoint ut_ini = BP_UNINIT, ut_start, ut_mid, ut_end, ut_endlast;
    BasePoint oxy, sxy;
    bigreal last_t, t;
    int is_right, linear, curved;
    int is_ccw_ini, is_ccw_start, is_ccw_mid, was_ccw;
    int closed = ss->first->prev!=NULL;

    if ( !c->remove_outer || !closed )
	left = chunkalloc(sizeof(SplineSet));
    if ( !c->remove_inner || !closed )
	right = chunkalloc(sizeof(SplineSet));

    for ( s=ss->first->next; s!=NULL && s!=first; s=s->to->next ) {
	if ( first==NULL )
	    first = s;

	if ( SplineLength(s)==0 )
	    continue;		// Can ignore zero length splines

	ut_start = SplineUTanVecAt(s, 0.0);
	linear = SplineIsLinearish(s);
	if ( linear ) {
	    ut_end = ut_start;
	    is_ccw_start = 0;
	} else {
	    ut_end = SplineUTanVecAt(s, 1.0);
	    is_ccw_start = SplineTurningCCWAt(s, 0.0);
	}
	if ( BP_IS_UNINIT(ut_ini) ) {
	    ut_ini = ut_start;
	    is_ccw_ini = is_ccw_start;
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
	    CalcNibOffset(c, ut_start, is_right, &no, -1);
	    oxy = BP_ADD(sxy, no.off[is_ccw_start]);

	    printf("is_right=%d nci=%d is_ccw:%d nc.x=%lf nc.y=%lf x=%lf, y=%lf, ut_start=%lf,%lf\n", is_right, no.nci[is_ccw_start], is_ccw_start, c->nibcorners[no.nci[is_ccw_start]].on_nib->me.x, c->nibcorners[no.nci[is_ccw_start]].on_nib->me.y, oxy.x, oxy.y, ut_start.x, ut_start.y);

	    // Create or verify initial spline point
	    if ( cur->first==NULL ) {
		cur->first = SplinePointCreate(oxy.x, oxy.y);
		cur->last = cur->first;
	    } else 
		HandleJoint(c, cur, sxy, oxy, &no, is_ccw_start, ut_endlast,
		            was_ccw, is_right);

	    // The path for this spline
	    if ( linear ) {
		sxy = SPLINEPVAL(s, 1.0);
		oxy = BP_ADD(sxy, no.off[is_ccw_start]);
		sp = SplinePointCreate(oxy.x, oxy.y);
		SplineMake3(cur->last, sp);
		cur->last = sp;
	    } else {
		t = 0.0;
		ut_mid = ut_start;
		is_ccw_mid = is_ccw_start;
		while ( t < 1.0 ) {
		    last_t = t;
		    t = SplineStrokeNextT(c, s, t, is_ccw_mid, &ut_mid, &curved,
		                          is_right, no.nci[is_ccw_mid]);
		    assert( t > last_t );
		    printf("nci: %d, ut_mid=%lf,%lf, t=%.15lf\n",
		           no.nci[is_ccw_mid], ut_mid.x, ut_mid.y, t);

		    if ( curved )
			sp = TraceAndFitSpline(c, s, last_t, t, cur->last,
			                       no.nci[is_ccw_mid], is_right);
		    else
			sp = AppendCubicSplinePortion(s, last_t, t, cur->last);

		    cur->last = sp;
		    sxy = SPLINEPVAL(s, t);
		    CalcNibOffset(c, ut_mid, is_right, &no, no.nci[is_ccw_mid]);
		    SplineStrokeAppendFixup(cur->last, sxy, &no);

		    // Deal with possible corner ambiguity
		    if (t < 1.0) {
			is_ccw_mid = SplineTurningCCWAt(s, t);
			oxy = BP_ADD(sxy, no.off[is_ccw_mid]);
			if ( !BPNEAR(cur->last->me, oxy) ) {
			    assert(BPNEAR(cur->last->me,
			           (BP_ADD(sxy, no.off[!is_ccw_mid]))));
			    sp = SplinePointCreate(oxy.x, oxy.y);
			    SplineMake3(cur->last, sp);
			    cur->last = sp;
			}
		    }
		}
	    }
	}
	ut_endlast = ut_end;
        was_ccw = SplineTurningCCWAt(s, 1.0);
    }
    if (    (left!=NULL && left->first==NULL) 
         || (right!=NULL && right->first==NULL) ) {
	// Presumably the path had only zero-length splines
	assert(    (left==NULL || left->first==NULL)
	        && (right==NULL || right->first==NULL) );
	chunkfree(left, sizeof(SplineSet));
	chunkfree(right, sizeof(SplineSet));
	return NULL;
    }
    if ( !closed ) {
	CalcNibOffset(c, ut_endlast, false, &no_last, -1);
	CalcNibOffset(c, ut_endlast, true, &no, -1);
	oxy = BP_ADD(ss->last->me, no_last.off[was_ccw]);
	if ( !BPNEAR(left->last->me, oxy) ) {
	    assert(BPNEAR(left->last->me,
	                  (BP_ADD(ss->last->me, no_last.off[!was_ccw]))));
	    sp = SplinePointCreate(oxy.x, oxy.y);
	    SplineMake3(left->last, sp);
	    left->last = sp;
	}
        sp = AddJoint(c, left->last, &no_last, was_ccw, &no, was_ccw, false);
	left->last = sp;
	SplineSetReverse(right);
	left->next = right;
	right = NULL;
	SplineSetJoin(left, true, FIXUP_MARGIN, &closed);
	if ( !closed )
	     LogError( _("Warning: Contour did not close as expected\n") );
	CalcNibOffset(c, ut_ini, true, &no_last, -1);
	CalcNibOffset(c, ut_ini, false, &no, -1);
        sp = AddJoint(c, left->last, &no, is_ccw_start, &no_last,
	              is_ccw_start, true);
	left->last = sp;
	SplineSetJoin(left, true, FIXUP_MARGIN, &closed);
	if ( !closed )
	     LogError( _("Warning: Contour did not close as expected\n") );
	else if ( c->rmov==srmov_contour )
	    left = SplineSetRemoveOverlap(NULL,left,over_remove);
    } else {
	// This can fail if the source contours are closed in a strange way
	if ( left!=NULL ) {
	    if ( !BPWITHIN(left->first->me, left->last->me,
	                   INTERSPLINE_MARGIN) ) {
		CalcNibOffset(c, ut_ini, false, &no, -1);
		HandleJoint(c, left, sxy, left->first->me, &no, is_ccw_ini,
		            ut_endlast, was_ccw, false);
	    }
            left = SplineSetJoin(left, true, INTERSPLINE_MARGIN, &closed);
	    if ( !closed )
		LogError( _("Warning: Left contour did not close\n") );
	    else if ( c->rmov==srmov_contour )
		left = SplineSetRemoveOverlap(NULL,left,over_remove);
	}
	if ( right!=NULL ) {
	    if ( !BPWITHIN(right->first->me, right->last->me,
	                   INTERSPLINE_MARGIN) ) {
		CalcNibOffset(c, ut_ini, true, &no, -1);
		HandleJoint(c, right, sxy, right->first->me, &no, is_ccw_ini,
		            ut_endlast, was_ccw, true);
	    }
            right = SplineSetJoin(right, true, INTERSPLINE_MARGIN, &closed);
	    if ( !closed )
		LogError( _("Warning: Right contour did not close\n") );
	    else if ( c->rmov!=srmov_none ) {
		SplineSetReverse(right);
		// Need to do this for either srmov_contour or srmov_layer
		right = SplineContourOuterCCWRemoveOverlap(right);
	    }
	    if ( left != NULL ) {
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

static SplinePointList *SinglePointStroke(SplinePoint *sp,
                                          struct strokecontext *c) {
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
	SplinePointListTransformExtended(ret, trans, tpt_AllPoints,
	                                 tpmask_dontTrimValues);
    }
    return( ret );
}

static SplineSet *SplineSet_Stroke(SplineSet *ss,struct strokecontext *c,
                                   int order2) {
    SplineSet *base = ss;
    SplineSet *ret;
    int copied = 0;
    int is_ccw = base->first->prev!=NULL && !SplinePointListIsClockwise(base);

    if ( base->first->next==NULL )
	ret = SinglePointStroke(base->first, c);
    else {
	if ( base->first->next->order2 ) {
	    base = SSPSApprox(ss);
	    copied = 1;
	}
	if ( is_ccw )
	    SplineSetReverse(base);
	ret = SplinesToContours(base, c);
    }
    if ( order2 )
	ret = SplineSetsConvertOrder(ret,order2);
    if ( copied )
	SplinePointListFree(base);
    else if ( is_ccw )
	SplineSetReverse(base);

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
    int max_pc;
    StrokeContext c;
    SplineSet *nibs, *nib, *bnext, *first, *last, *cur;
    bigreal sn = 0.0, co = 1.0; // zero degree defaults
    DBounds b;
    real trans[6];
    static struct simplifyinfo smpl = { sf_normal, 0.75, 0.2, 0, 0, 0, 0 };

    if ( si->stroke_type==si_centerline )
	IError("centerline not handled");

    memset(&c,0,sizeof(c));
    c.join = si->join;
    c.cap  = si->cap;
    c.miterlimit = /* -cos(theta) */ -.98 /* ~11 degrees, PS miterlimit=10 */;
    c.remove_inner = si->removeinternal;
    c.remove_outer = si->removeexternal;
    c.leave_users_center = si->leave_users_center;
    c.extrema = !si->noextrema;
    c.simplify = !si->nosimplify;
    c.rmov = si->rmov;
    //c.scaled_or_rotated = si->factor!=NULL;

    c.leave_users_center = 1;

    if ( si->penangle!=0 ) {
	sn = sin(si->penangle);
	co = cos(si->penangle);
    }
    trans[0] = co;
    trans[1] = sn;
    trans[2] = -sn;
    trans[3] = co;
    trans[4] = trans[5] = 0;

    if ( si->stroke_type==si_std || si->stroke_type==si_caligraphic ) {
	c.pentype = si->stroke_type==si_std ? pt_circle : pt_square;
	max_pc = 4;
	nibs = UnitShape(si->stroke_type==si_std ? 0 : -4);
	trans[0] *= si->radius;
	trans[1] *= si->radius;
	trans[2] *= si->minorradius;
	trans[3] *= si->minorradius;
    } else {
	c.pentype = pt_convex;
	max_pc = 20; // a guess
	nibs = SplinePointListCopy(si->nib);
	if ( !c.leave_users_center ) {
	    SplineSetQuickBounds(nibs,&b);
	    trans[4] = -(b.minx+b.maxx)/2;
	    trans[5] = -(b.miny+b.maxy)/2;
	}
    }
    SplinePointListTransformExtended(nibs,trans,tpt_AllPoints,
	                             tpmask_dontTrimValues);

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
    if ( c.rmov==srmov_layer )
	first=SplineSetRemoveOverlap(NULL,first,over_remove);
    if ( c.extrema )
	SplineCharAddExtrema(NULL,first,ae_all,1000);
    if ( c.simplify )
	first = SplineCharSimplify(NULL,first,&smpl);
    else
	SPLCategorizePoints(first);

    free(c.nibcorners);
    SplinePointListFree(nibs);
    nibs = NULL;

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

void SplineStrokeTests() {
    int i, j, k, a, b;
    bigreal r, x, y, z;
    BasePoint ut[361], utmax = { -1, 0 }, utmin = UTMIN;

    printf("Spline Stroke Tests\n");

    for (i = 0; i<=360; i++) {
	r = ((180-i) * PI / 180);
	ut[i] = UTanVectorize(cos(r), sin(r));
//	printf("i:%d, i.x:%lf, i.y:%lf\n", i, ut[i].x, ut[i].y);
    }

    for (i=0; i<360; ++i) {
	for (j = 0; j<360; ++j) {
	    if ( i==j )
		continue;
	    if ( i<j )
		assert(    UTanVecGreater(ut[i], ut[j])
		       && !UTanVecGreater(ut[j], ut[i]));
	    else
		assert(    UTanVecGreater(ut[j], ut[i])
		       	&& !UTanVecGreater(ut[i], ut[j]));
	    assert(!BPNEAR(ut[i], ut[j]));
	    x = atan2(ut[i].x, ut[i].y);
	    y = atan2(ut[j].x, ut[j].y);
	    z = x-y;
	    z = atan2(sin(z), cos(z));
	    if ( RealNear(z, PI) || RealNear(z, -PI) )
		continue;
	    assert( (z>0) == JointBendsCW(ut[i], ut[j]) );
	}
	assert( i==0 || UTanVecGreater(utmax, ut[i]) );
	assert( UTanVecGreater(ut[i], utmin) );
    }

    for (i=0; i<360; ++i)
	for (j = 0; j<360; ++j)
	    for (k = 0; k<360; ++k) {
		a = UTanVecsSequent(ut[i], ut[j], ut[k], false);
		b = UTanVecsSequent(ut[i], ut[j], ut[k], true);
		if ( i==j )
		    assert(a && b);
		else if ( j==k || i==k )
		    assert(!a && !b);
		else {
		    x = atan2(ut[i].x, ut[i].y);
		    y = atan2(ut[j].x, ut[j].y);
		    z = atan2(ut[k].x, ut[k].y);
		    y = fmod(4*PI+y-x,2*PI);
		    z = fmod(4*PI+z-x,2*PI);
		    if ( y<z )
			assert( a && !b );
		    else
			assert( !a && b );
		}
	    }
}
