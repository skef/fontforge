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
#include "splinefit.h"
#include "splinefont.h"
#include "splineorder2.h"
#include "splineoverlap.h"
#include "splineutil.h"
#include "splineutil2.h"
#include "utanvec.h"

#include <math.h>
#include <assert.h>

#define CIRCOFF 0.551915
#define PI      3.1415926535897932

// Rounding makes even slope-continuous splines only approximately so. 
#define INTERSPLINE_MARGIN (1e-1)
#define INTRASPLINE_MARGIN (1e-8)
#define FIXUP_MARGIN (1e-2)
#define CUSPD_MARGIN (1e-5)
// About .25 degrees
#define COS_MARGIN (1e-5)
#define MIN_ACCURACY (1e-5)

#define NORMANGLE(a) ((a)>PI?(a)-2*PI:(a)<-PI?(a)+2*PI:(a))
#define BPNEAR(bp1, bp2) BPWITHIN(bp1, bp2, INTRASPLINE_MARGIN)
#define BP_DIST(bp1, bp2) sqrt(pow((bp1).x-(bp2).x,2)+pow((bp1).y-(bp2).y,2))

enum pentype { pt_circle, pt_square, pt_convex };

/* c->nibcorners is a structure that models the "corners" of the current nib
 * treated as if all splines were linear, and therefore as a convex polygon.
 */

#define NC_IN_IDX 0
#define NC_OUT_IDX 1

#define N_NEXTI(n, i) (((i)+1)%(n))
#define N_PREVI(n, i) (((n)+(i)-1)%(n))
#define NC_NEXTI(c, i) N_NEXTI((c)->n, i)
#define NC_PREVI(c, i) N_PREVI((c)->n, i)

typedef struct nibcorner {
    SplinePoint *on_nib;
    BasePoint utv[2];    // Unit tangents entering and leaving corner
    unsigned int linear: 1;
} NibCorner;

// XXX Add optional error info char *
typedef struct strokecontext {
    enum pentype pentype;
    enum linejoin join;
    enum linecap cap;
    enum stroke_rmov rmov;
    bigreal joinlimit;
    bigreal extendcap;
    bigreal acctarget;
    SplineSet *nib;
    int n;
    NibCorner *nibcorners;
    BasePoint pseudo_origin;
    unsigned int jlrelative: 1;
    unsigned int ecrelative: 1;
    unsigned int remove_inner: 1;
    unsigned int remove_outer: 1;
    unsigned int leave_users_center: 1;	// Don't center the nib when stroking
    unsigned int simplify:1;
    unsigned int extrema:1;
    unsigned int contour_was_ccw:1;
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

StrokeInfo *InitializeStrokeInfo(StrokeInfo *sip) {
    if ( sip==NULL )
	sip = malloc(sizeof(StrokeInfo));

    *sip = (StrokeInfo) {
        25.0,         // radius
        lj_nib,       // line_join
        lc_nib,       // line_cap
        si_round,     // pen type
        srmov_layer,  // remove overlap style
        false,        // removeinternal
        false,        // removeexternal
	true,        // simplify
	true,        // extrema
        false,        // leave users center
        true,        // jlrelative
        true,        // ecrelative
        0.0,          // pen angle
        0.0,          // minor radius (0: major)
	0.0,          // extend cap
	20.0,         // join limit
	0.25,         // accuracy target
        NULL,         // nib (SplineSet for si_nib)
        0.0,          // freehand: radius 2
        0.0,          // freehand: pressure 1
        0.0,          // freehand: pressure 2
        NULL,         // data
        NULL          // factor
    };

    return sip;
}

void SITranslatePSArgs(StrokeInfo *sip, enum linejoin lj, enum linecap lc) {
    sip->extendcap = 0;
    switch (lc) {
	case lc_round:
	    sip->cap = lc_nib;
	    break;
	case lc_square:
	    sip->cap = lc_butt;
	    sip->extendcap = 1;
	    sip->ecrelative = true;
	    break;
	default:
	    sip->cap = lc;
    }
    switch (lj) {
	case lc_round:
	    sip->join = lj_nib;
	    break;
	default:
	    sip->join = lj;
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

int ConvexNibID(char *tok) {
    if ( tok!=NULL ) {
	if ( strcmp(tok, "default") )
	    return 1;
	else if ( strcmp(tok, "ui") )
	    return -1;
	else if ( strcmp(tok, "freehand") )
	    return -2;
	else {
	    return 0;
	}
    } else
	return 1;
}

static SplineSet *convex_nibs[2];

int StrokeSetConvexPlain(SplineSet *ss, int toknum) {
    if ( toknum==1 ) {
	if ( convex_nibs[toknum]!=NULL )
	    SplinePointListFree(convex_nibs[toknum]);
	convex_nibs[toknum] = ss;
	return true;
    }

    return false;
}

SplineSet *StrokeGetConvexPlain(int toknum) {
    if ( toknum==1 )
	return convex_nibs[toknum];
    return NULL;
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
 */
// XXX Add outer CP check against next angle
enum ShapeType NibIsValid(SplineSet *ss) {
    Spline *s;
    int n = 1;
    bigreal d, anglesum = 0, angle, last_angle;
    bigreal pcp_angle, ncp_angle, last_pcp_angle;

    if ( ss->first==NULL )
	return Shape_TooFewPoints;
    if ( ss->first->prev==NULL )
	return Shape_NotClosed;
    if ( !SplinePointListIsClockwise(ss) )
	return Shape_CCW;

    s = ss->first->prev;
    last_angle = atan2(s->to->me.y - s->from->me.y,
                       s->to->me.x - s->from->me.x);
    if ( ! SplineIsLinear(s) )
	last_pcp_angle = atan2(s->to->me.y - s->to->prevcp.y,
	                       s->to->me.x - s->to->prevcp.x);
    else
	last_pcp_angle = last_angle;

    s = ss->first->next;
    SplinePointListSelect(ss, false);
    while ( true ) {
	s->from->selected = true;
	if ( BPWITHIN(s->from->me, s->to->me,1e-2) )
	    return Shape_TinySpline;
	angle = atan2(s->to->me.y - s->from->me.y, s->to->me.x - s->from->me.x);
	if ( RealWithin(angle, last_angle, 1e-4) )
	    return Shape_PointOnEdge;
	d = last_angle-angle;
	d = NORMANGLE(d);
	if ( d<0 )
	    return Shape_CCWTurn;
	anglesum += d;
	if ( last_pcp_angle != last_angle ) {
	    d = last_pcp_angle-angle;
	    d = NORMANGLE(d);
	    if ( d<0 )
		return Shape_BadPrevCP;
	}
	if ( !SplineIsLinear(s) ) {
	    if ( s->from->nonextcp || BPNEAR(s->from->nextcp, s->from->me) )
		return Shape_HalfLinear;
	    ncp_angle = atan2(s->from->nextcp.y - s->from->me.y,
	                      s->from->nextcp.x - s->from->me.x);
	    d = ncp_angle-last_pcp_angle; // Outer limit
	    d = NORMANGLE(d);
	    if ( d > 0.0 )
		return Shape_BadNextCP;
	    d = angle-ncp_angle; // Inner limit
	    d = NORMANGLE(d);
	    if ( d > 0.0 )
		return Shape_BadNextCP;
	    printf("nextcp angle: %lf, dout: %lf, din:%lf\n", ncp_angle, d, NORMANGLE(angle-ncp_angle));
	    s->from->selected = false;
	    s->to->selected = true;
	    if ( s->to->noprevcp || BPNEAR(s->to->prevcp, s->to->me) ) 
		return Shape_HalfLinear;
	    pcp_angle = atan2(s->to->me.y - s->to->prevcp.y,
	                      s->to->me.x - s->to->prevcp.x);
	    d = pcp_angle-angle; // Inner limit
	    d = NORMANGLE(d);
	    if ( d > 0.0 )
		return Shape_BadPrevCP;
	    printf("prevcp angle: %lf, dout: %lf\n", pcp_angle, d);
	    s->from->selected = false;
	} else
	    pcp_angle = angle;
	s->from->selected = s->to->selected = false;
	s=s->to->next;
	last_angle = angle;
	last_pcp_angle = pcp_angle;
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

static void BuildNibCorners(NibCorner **ncp, SplineSet *nib, int *maxp,
                            int *n) {
    int i, max_utan_index = -1;
    BasePoint max_utanangle = UTMIN;
    SplinePoint *sp;
    NibCorner *nc = *ncp, *tpc;

    if ( nc==NULL )
	nc = calloc(*maxp,sizeof(NibCorner));

    for ( sp=nib->first, i=0; ; ) {
	if ( i==*maxp ) { // We guessed wrong
	    *maxp *= 2;
	    tpc = calloc(*maxp, sizeof(NibCorner));
	    memcpy(tpc, nc, (i-1)*sizeof(NibCorner));
	    free(nc);
	    nc = tpc;
	}
	nc[i].on_nib = sp;
	nc[i].linear = SplineIsLinear(sp->next);
	nc[i].utv[NC_IN_IDX] = SplineUTanVecAt(sp->prev, 1.0);
	nc[i].utv[NC_OUT_IDX] = SplineUTanVecAt(sp->next, 0.0);
	if (    UTanVecGreater(nc[i].utv[NC_IN_IDX], max_utanangle)
	     || BPNEAR(nc[i].utv[NC_IN_IDX], max_utanangle) ) {
	    max_utan_index = i;
	    max_utanangle = nc[i].utv[NC_IN_IDX];
	}
	++i;
	sp=sp->next->to;
	if ( sp==nib->first )
	    break;
    }
    *n = i;

    // Put in order of decreasing utanvec[NC_IN_IDX]
    if (max_utan_index != 0) {
	tpc = malloc(i*sizeof(NibCorner));
	memcpy(tpc, nc+max_utan_index, (i-max_utan_index)*sizeof(NibCorner));
	memcpy(tpc+(i-max_utan_index), nc, max_utan_index*sizeof(NibCorner));
	free(nc);
	nc = tpc;
    }

    *ncp = nc;
}

// The index into c->nibcorners associated with unit tangent ut
static int _IndexForUTanVec(NibCorner *nc, int n, BasePoint ut, int nci_hint) {
    int nci;

    assert( nci_hint == -1 || (nci_hint >= 0 && nci_hint < n ) );

    if (   nci_hint != -1
        && UTanVecsSequent(nc[nci_hint].utv[NC_IN_IDX], ut,
	                   nc[N_NEXTI(n, nci_hint)].utv[NC_IN_IDX],
                           false) )
	nci = nci_hint;
    else // XXX Replace with binary search
	for (nci = 0; nci<n; ++nci)
	    if ( UTanVecsSequent(nc[nci].utv[NC_IN_IDX], ut,
	                         nc[N_NEXTI(n, nci)].utv[NC_IN_IDX],
	                         false) )
		break;
    return nci;
}

static int IndexForUTanVec(StrokeContext *c, BasePoint ut, int nci_hint) {
    return _IndexForUTanVec(c->nibcorners, c->n, ut, nci_hint);
}

/* The offset from the stroked curve associated with unit tangent ut
 * (or it's reverse).
 */
static NibOffset *_CalcNibOffset(NibCorner *nc, int n, BasePoint ut,
                                 int reverse, NibOffset *no, int nci_hint) {
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

    nci = no->nci[0] = no->nci[1] = _IndexForUTanVec(nc, n, ut, nci_hint);
    ncpi = N_PREVI(n, nci);
    ncni = N_NEXTI(n, nci);

    // The two cases where one of two points might draw the same angle
    // and therefore the only cases where the array values differ
    if (   nc[nci].linear
	// "Capsule" case
        && BPNEAR(ut, nc[ncni].utv[NC_IN_IDX])
        && BPNEAR(ut, nc[nci].utv[NC_IN_IDX]) ) {
	no->nt = 0.0;
	no->off[NIBOFF_CCW_IDX] = nc[nci].on_nib->me;
	no->off[NIBOFF_CW_IDX] = nc[ncni].on_nib->me;
	no->nci[NIBOFF_CW_IDX] = ncni;
	no->at_line = true;
	no->curve = false;
    } else if ( nc[ncpi].linear && BPNEAR(ut, nc[ncpi].utv[NC_OUT_IDX]) ) {
	// Other lines
	no->nt = 0.0;
	no->off[NIBOFF_CCW_IDX] = nc[ncpi].on_nib->me;
	no->nci[NIBOFF_CCW_IDX] = ncpi;
	no->off[NIBOFF_CW_IDX] = nc[nci].on_nib->me;
	no->at_line = true;
	no->curve = false;
    // Whether ut is (effectively) between IN and OUT, in 
    // which case the point draws the offset curve
    } else if (   UTanVecsSequent(nc[nci].utv[NC_IN_IDX], ut,
                                  nc[nci].utv[NC_OUT_IDX], false) 
               || BPNEAR(ut, nc[nci].utv[NC_OUT_IDX] ) ) {
	no->nt = 0.0;
	no->off[0] = no->off[1] = nc[nci].on_nib->me;
	no->curve = false;
    } else {
	// ut is on a spline curve clockwise from point nci
	assert( UTanVecsSequent(nc[nci].utv[NC_OUT_IDX], ut,
	                        nc[ncni].utv[NC_IN_IDX], false) );
	// Nib splines are locally convex and therefore have t value per slope
	ns = nc[nci].on_nib->next;
	no->nt = SplineSolveForUTanVec(ns, ut, 0.0);
	no->off[0] = no->off[1] = SPLINEPVAL(ns, no->nt);
	no->curve = true;
    }
    return no;
}

static NibOffset *CalcNibOffset(StrokeContext *c, BasePoint ut, int reverse,
                                NibOffset *no, int nci_hint) {
    return _CalcNibOffset(c->nibcorners, c->n, ut, reverse, no, nci_hint);
}

/* This is a convenience function that also concisely demonstrates
 * the calculation of an offseted point
 */
BasePoint SplineOffsetAt(StrokeContext *c, Spline *s, bigreal t, int is_right) {
    int is_ccw;
    BasePoint xy, ut;
    NibOffset no;

    xy = SPLINEPVAL(s, t);
    is_ccw = SplineTurningCCWAt(s, t);

    ut = SplineUTanVecAt(s, t);
    CalcNibOffset(c, ut, is_right, &no, -1);

    return BP_ADD(xy, no.off[is_ccw]);
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
    if ( RealWithin(t_start, t_end, 1e-4) )
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
 * after point dst_start. "backward" (and therefore forward) is relative
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
	assert( s!=NULL && s!=s_start ); // XXX turn into runtime warning?
        dst_start = AppendCubicSplinePortion(s, backward ? 1 : 0,
	                                     backward ? 0 : 1, dst_start);
    }
    dst_start = AppendCubicSplinePortion(s, backward ? 1 : 0, t_end, dst_start);
    return dst_start;
}

static SplinePoint *AddNibPortion(NibCorner *nc, SplinePoint *start_p,
                                  NibOffset *nos, int is_ccw_s, NibOffset *noe,
			          int is_ccw_e, int bk) {
    SplinePoint *sp, *start_np, *end_np;

    start_np = nc[nos->nci[is_ccw_s]].on_nib;
    end_np = nc[noe->nci[is_ccw_e]].on_nib;
    sp = AppendCubicSplineSetPortion(start_np->next, nos->nt, end_np->next,
                                     noe->nt, start_p, bk);
    return sp;
}

static void SSAppendSemiCircle(SplineSet *cur, bigreal radius, BasePoint ut,
                               int bk) {
    SplineSet *circ;
    real trans[6];
    BasePoint cxy = cur->last->me;
    SplinePoint *sp;
    NibOffset nos, noe;
    NibCorner *nc = NULL;
    int n, mn=4;

    circ = UnitShape(0);
    memset(&trans, 0, sizeof(trans));
    trans[0] = trans[3] = radius;
    SplinePointListTransformExtended(circ, trans, tpt_AllPoints,
	                             tpmask_dontTrimValues);
    BuildNibCorners(&nc, circ, &mn, &n);
    assert( n==4 && nc!=NULL );

    _CalcNibOffset(nc, n, ut, false, &nos, -1);
    _CalcNibOffset(nc, n, ut, true, &noe, -1);
    sp = AddNibPortion(nc, cur->last, &nos, false, &noe, false, bk);
    cur->last = sp;
    SplinePointListFree(circ);
    free(nc);
}

/* Put the new endpoint exactly where the NibOffset calculation says it
 * should be to avoid cumulative append errors.
 */
static void SplineStrokeAppendFixup(SplinePoint *endp, BasePoint sxy,
                                    NibOffset *noe, int is_ccw) {
    int i = 0;
    bigreal mg;
    BasePoint oxy, dxy;

    oxy = BP_ADD(sxy, noe->off[is_ccw]);
    dxy = BP_ADD(oxy, BP_REV(endp->me));
    mg = fmax(fabs(dxy.x), fabs(dxy.y));

    assert( mg < 1 );
    if ( mg > FIXUP_MARGIN ) {
	LogError(_("Warning: Coordinate diff %lf greater than margin %lf\n"),
	         mg, FIXUP_MARGIN);
    }

    endp->nextcp = BP_ADD(endp->nextcp, dxy);
    endp->me = oxy;
}

/* Note that is_ccw relates to the relevant NibOffset entries, and not
 * necessarily the direction of the curve at t.
 */
static int OffsetOnCuspAt(StrokeContext *c, Spline *s, bigreal t,
                          NibOffset *nop, int is_right, int is_ccw)
{
    bigreal cs, cn;
    NibOffset no;

    cs = SplineCurvature(s, t);

    // Cusps aren't a problem when increasing the radius of the curve
    if ( (cs>0) == !!is_right )
	return false;

    if ( nop==NULL ) {
	nop = &no;
	CalcNibOffset(c, SplineUTanVecAt(s, t), is_right, nop, -1);
    }
    cn = SplineCurvature(c->nibcorners[nop->nci[is_ccw]].on_nib->next, nop->nt);

    // printf("t: %lf, SplineCurvature: %lf and %lf, right: %d, ccw: %d\n", t, cs, cn, is_right, is_ccw);
    return RealWithin(cn, 0, 1e-6) ? false : is_right ? (cn >= cs) : (cn >= -cs);
}

static bigreal SplineFindCuspSing(StrokeContext *c, Spline *s, bigreal t_start,
                                  bigreal t_end, int is_right, int is_ccw,
				  bigreal margin, int on_cusp)
{
    bigreal t_mid;
    int cusp_mid;
    assert( t_start >= 0.0 && t_start <= 1.0 );
    assert( t_end >= 0.0 && t_end <= 1.0 );
    assert( t_start < t_end );
    assert( OffsetOnCuspAt(c, s, t_start, NULL, is_right, is_ccw) == on_cusp );
    assert( OffsetOnCuspAt(c, s, t_end, NULL, is_right, is_ccw) != on_cusp );
    while ( t_end-t_start > margin ) {
	t_mid = (t_start+t_end)/2;
	cusp_mid = OffsetOnCuspAt(c, s, t_mid, NULL, is_right, is_ccw);
	if ( cusp_mid==on_cusp )
	    t_start = t_mid;
	else
	    t_end = t_mid;
    }
    return t_end;
}

typedef struct stroketraceinfo {
    StrokeContext *c;
    Spline *s;
    bigreal cusp_trans;
    int nci_hint;
    int num_points;
    unsigned int is_right: 1;
    unsigned int starts_on_cusp: 1;
    unsigned int first_pass: 1;
    unsigned int found_trans: 1;
} StrokeTraceInfo;

int GenStrokeTracePoints(void *vinfo, bigreal t_start, bigreal t_end,
                         FitPoint **fpp) {
    StrokeTraceInfo *stip = (StrokeTraceInfo *)vinfo;
    int i, nib_ccw, on_cusp;
    NibOffset no;
    FitPoint *fp, tmp_fp;
    bigreal nidiff, t, csd;
    BasePoint xy, rel_ut;

    fp = calloc(stip->num_points, sizeof(FitPoint));
    nidiff = (t_end - t_start) / (stip->num_points-1);

    nib_ccw = SplineTurningCCWAt(stip->s, t_start);
    no.nci[0] = no.nci[1] = stip->nci_hint;
    for ( i=0, t=t_start; i<stip->num_points; ++i, t+= nidiff ) {
	if ( i==(stip->num_points-1) ) {
	    nib_ccw = !nib_ccw; // Stop at the closer corner
	    t = t_end; // side-step nidiff rounding errors
	}
	xy = SPLINEPVAL(stip->s, t);
	fp[i].ut = SplineUTanVecAt(stip->s, t);
	CalcNibOffset(stip->c, fp[i].ut, stip->is_right, &no, no.nci[nib_ccw]);
	fp[i].p = BP_ADD(xy, no.off[nib_ccw]);
	fp[i].t = t;
	if ( stip->first_pass ) {
	    on_cusp = OffsetOnCuspAt(stip->c, stip->s, t, &no,
	                             stip->is_right, nib_ccw);
	    if ( on_cusp!=stip->starts_on_cusp ) {
		stip->found_trans = true;
		stip->cusp_trans = SplineFindCuspSing(stip->c, stip->s, t-nidiff,
		                                      t, stip->is_right, 
		                                      nib_ccw, CUSPD_MARGIN,
						      stip->starts_on_cusp);
		free(fp);
		printf("Found %s cusp transition at %lf\n", stip->starts_on_cusp ? "1 to 0" : "0 to 1", stip->cusp_trans);
		return 0;
	    }
	} else {
#ifndef NODEBUG
	;
#else 
	;
#endif // NODEBUG
	}
	if ( stip->starts_on_cusp )
	    fp[i].ut = BP_REV(fp[i].ut);
    }
    *fpp = fp;
    stip->first_pass = false;
    return stip->num_points;
}

#define TRACE_CUSPS false
SplinePoint *TraceAndFitSpline(StrokeContext *c, Spline *s, bigreal t_start,
                               bigreal t_end, SplinePoint *start,
			       int nci_hint, int is_right, int on_cusp) {
    SplinePoint *sp = NULL;
    StrokeTraceInfo sti;
    FitPoint *fpp;
    BasePoint xy;

    sti.c = c;
    sti.s = s;
    sti.nci_hint = nci_hint;
    sti.num_points = 10;
    sti.is_right = is_right;
    sti.first_pass = true;
    sti.starts_on_cusp = on_cusp;
    sti.found_trans = false;

    if ( on_cusp && !TRACE_CUSPS )
	GenStrokeTracePoints((void *)&sti, t_start, t_end, &fpp);
    else 
	sp = ApproximateSplineSetFromGen(start, NULL, t_start, t_end, 
	                                 c->acctarget, false,
	                                 &GenStrokeTracePoints, (void *) &sti,
	                                 false);
    if ( !sti.found_trans ) {
	if ( sp==NULL ) {
	    assert( on_cusp && !TRACE_CUSPS );
	    xy = SplineOffsetAt(c, s, t_end, is_right);
	    sp = SplinePointCreate(xy.x, xy.y);
	    SplineMake3(start, sp);
	}
   	return sp;
    }

    assert ( sti.found_trans && sp==NULL );
    assert ( sti.cusp_trans >= t_start && sti.cusp_trans <= t_end );
    if ( on_cusp && !TRACE_CUSPS ) {
	xy = SplineOffsetAt(c, s, sti.cusp_trans, is_right);
	sp = SplinePointCreate(xy.x, xy.y);
	SplineMake3(start, sp);
    } else {
	sti.first_pass = false;
	sp = ApproximateSplineSetFromGen(start, NULL, t_start, sti.cusp_trans, 
	                                 c->acctarget, false,
	                                 &GenStrokeTracePoints, (void *) &sti,
	                                 false);
    }
    assert( sp!=NULL );
    sp->pointtype = pt_corner;
    on_cusp = !on_cusp;
    if ( RealWithin(sti.cusp_trans, t_end, CUSPD_MARGIN) )
	return sp;
    return TraceAndFitSpline(c, s, sti.cusp_trans, t_end, sp, nci_hint,
                             is_right, on_cusp);
}
#undef TRACE_CUSPS

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
		if ( BPNEAR(c->nibcorners[nci].utv[NC_IN_IDX],
		            c->nibcorners[ncpi].utv[NC_IN_IDX]) ) {
		    assert(c->nibcorners[ncpi].linear);
		    *curved = true;
		    inout = NC_OUT_IDX;
		    nci = NC_PREVI(c, ncpi);
		} else {
		    assert(c->nibcorners[ncpi].linear);
		    inout = NC_IN_IDX;
		    *curved = false;
		    nci = ncpi;
		}
	    } else {
		inout = NC_OUT_IDX;
		*curved = true;
		nci = ncpi;
	    }
	} else { // CW
	    if ( BPNEAR(c->nibcorners[nci].utv[NC_IN_IDX],
	                c->nibcorners[nci].utv[NC_OUT_IDX]) ) {
		if ( BPNEAR(c->nibcorners[nci].utv[NC_IN_IDX],
		            c->nibcorners[ncni].utv[NC_IN_IDX]) ) {
		    assert(c->nibcorners[nci].linear);
		    inout = NC_OUT_IDX;
		    *curved = true;
		} else {
		    inout = NC_IN_IDX;
		    *curved = !c->nibcorners[nci].linear;
		}
		nci = ncni;
	    } else {
		inout = NC_OUT_IDX;
		*curved = false;
	    }
	}
    } else if ( BPNEAR(ut, c->nibcorners[nci].utv[NC_OUT_IDX]) ) {
	assert( ! BPNEAR(c->nibcorners[nci].utv[NC_IN_IDX],
	                 c->nibcorners[nci].utv[NC_OUT_IDX]) );
	if ( is_ccw ) {
	    inout = NC_IN_IDX;
	    *curved = false;
	} else { // CW
	    if ( BPNEAR(c->nibcorners[nci].utv[NC_OUT_IDX],
	                c->nibcorners[ncni].utv[NC_IN_IDX]) ) {
		assert(c->nibcorners[nci].linear);
		inout = NC_OUT_IDX;
		*curved = !c->nibcorners[ncni].linear;
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

    // If there is an inflection point before next_t the spline will start
    // curving in the opposite direction, so stop and trace the next section
    // separately. An alternative would be find the next stop angle in the
    // other direction by returning:
    //
    // return SplineStrokeNextT(c, s, poi[i], !is_ccw, cur_ut, curved,
    //                          reverse, nci_hint);
    //
    // This would have one main advantage, which is that the offset inflection
    // points will not exactly correspond to the position on s. Therefore
    // stopping at the s inflection leads to an offset point frustratingly close
    // to but not at the inflection point. The disadvantage is that we would 
    // need to track how many times the direction changes to feed the right
    // value to SplineStrokeAppendFixup, leading to more code complexity. 
    if ( icnt = Spline2DFindPointsOfInflection(s, poi) ) {
	assert ( icnt < 2 || poi[0] <= poi[1] );
	for ( i=0; i<2; ++i )
	    if (    poi[i] > cur_t
	         && !RealNear(poi[i], cur_t)
	         && (next_t==-1 || poi[i] < next_t) ) {
		next_t = poi[i];
		next_ut = SplineUTanVecAt(s, next_t);
		break;
	    }
    }
    if ( RealWithin(next_t, 1.0, 1e-4) ) // XXX distance-based?
	next_t = 1.0;

    if ( next_t==-1 ) {
	next_t = 1.0;
	*cur_ut = SplineUTanVecAt(s, next_t);
    } else
	*cur_ut = next_ut;

    *curved = next_curved;
    return next_t;
}

static void SplineSetLineTo(SplineSet *cur, BasePoint xy) {
    SplinePoint *sp = SplinePointCreate(xy.x, xy.y);
    SplineMake3(cur->last, sp);
    cur->last = sp;
}

static void HandleFlat(SplineSet *cur, BasePoint sxy, NibOffset *noi,
                       int is_ccw) {
    BasePoint oxy;
    SplinePoint *sp;

    oxy = BP_ADD(sxy, noi->off[is_ccw]);
    if ( !BPNEAR(cur->last->me, oxy) ) {
	assert( BPNEAR(cur->last->me, (BP_ADD(sxy, noi->off[!is_ccw]))) );
	SplineSetLineTo(cur, oxy);
    }
}

static bigreal FalseStrokeWidth(StrokeContext *c, BasePoint ut) {
    SplineSet *tss;
    DBounds b;
    real trans[6];

    assert( RealNear(BP_LENGTHSQ(ut),1) );
    // Rotate nib and grab height as span
    trans[0] = trans[3] = ut.x;
    trans[1] = ut.y;
    trans[2] = -ut.y;
    trans[4] = trans[5] = 0;
    tss = SplinePointListCopy(c->nib);
    SplinePointListTransformExtended(tss, trans, tpt_AllPoints,
	                             tpmask_dontTrimValues);
    SplineSetFindBounds(tss,&b);
    SplinePointListFree(tss);
    // printf("maxy: %lf, miny: %lf\n", b.maxy, b.miny);
    return b.maxy - b.miny;
}

static bigreal CalcJoinLength(bigreal fsw, BasePoint ut1, BasePoint ut2) {
    bigreal costheta = BP_DOT(ut1, ut2);

    // Angle of interest is pi - theta, so formula is fsw/sin((pi - theta)/2)
    //   = sin(pi/2 - theta/2) = sin(pi/2)cos(theta/2) - cos(pi/2)sin(theta/2)
    //   = 1 * cos(theta/2) - 0 * sin(theta/2) = cos(theta/2)
    //   = sqrt((1 + costheta)/2)
    return fsw / sqrt((1 + costheta)/2);
}

static bigreal CalcJoinLimit(StrokeContext *c, bigreal fsw) {
    if (c->joinlimit <= 0)
	return DBL_MAX;

    if (!c->jlrelative)
	return c->joinlimit;

    return c->joinlimit * fsw / 2;
}

static bigreal CalcCapExtend(StrokeContext *c, bigreal cw) {
    if (c->extendcap <= 0)
	return 0;

    if (!c->ecrelative)
	return c->extendcap;

    return c->extendcap * cw / 2;
}

static int LineSameSide(BasePoint l1, BasePoint l2, BasePoint p, BasePoint r) {
    bigreal t1, t2;

    t1 = (l2.x-l1.x)*(r.y-l1.y) - (r.x-l1.x)*(l2.y-l1.y);
    t2 = (l2.x-l1.x)*(p.y-l1.y) - (p.x-l1.x)*(l2.y-l1.y);

    return    !RealWithin(t1, 0, 1e-4) && !RealWithin(t2, 0, 1e-4)
           && (t1>0)==(t2>0);
}

static bigreal LineDist(BasePoint l1, BasePoint l2, BasePoint p) {
    assert (l1.x!=l2.x || l1.y!=l2.y);
    return   abs((l2.y-l1.y)*p.x - (l2.x-l1.x)*p.y + l2.x*l1.y -l2.y*l1.x)
           / BP_DIST(l1, l2);
}

/*
static void CalcDualBasis(BasePoint v1, BasePoint v2, BasePoint *q1,
                          BasePoint *q2) {
    if ( v1.x==0 ) {
	assert( v1.y!=0 && v2.x!=0 );
	q2->y = 0;
	q2->x = 1/v2.x;
	q1->y = 1/v1.y;
	q1->x = -v2.y / ( v1.y * v2.x );
    } else if ( v1.y==0 ) { 
	assert( v2.y!=0 );
	q2->x = 0;
	q2->y = 1/v2.y;
	q1->x = 1/v1.x;
	q1->y = -v2.x / ( v1.x * v2.y );
    } else if ( v2.x==0 ) {
	assert( v2.y!=0 );
	q1->y = 0;
	q1->x = 1/v1.x;
	q2->y = 1/v2.y;
	q2->x = -v1.y / ( v2.y * v1.x );
    } else if ( v2.y==0 ) { 
	q1->x = 0;
	q1->y = 1/v1.y;
	q2->x = 1/v2.x;
	q2->y = -v1.x / ( v2.x * v1.y );
    } else {
	q1->y = 1 / (v1.y - (v1.x * v2.y / v2.x));
	q2->y = 1 / (v2.y - (v2.x * v1.y / v1.x));
	q1->x = -v2.y * q1->y / v2.x;
	q2->x = -v1.y * q2->y / v1.x;
    }
    assert(    RealNear(BP_DOT(v1,*q1), 1) && RealNear(BP_DOT(v2,*q2), 1)
            && RealNear(BP_DOT(v1,*q2), 0) && RealNear(BP_DOT(v2,*q1), 0) );
}
*/

static void CalcExtend(BasePoint refp, BasePoint ut, BasePoint op1,
                       BasePoint op2, bigreal min,
                       BasePoint *np1, BasePoint *np2) {
    bigreal extra=0, tmp;
    int intersects;
    BasePoint cop1 = BP_ADD(op1, ut), cop2 = BP_ADD(op2, ut);
    BasePoint clip1 = BP_ADD(refp, BP_SCALE(ut, min));
    BasePoint clip2 = BP_ADD(clip1, UT_90CW(ut));

    if ( !LineSameSide(clip1, clip2, op1, refp) )
	extra = LineDist(clip1, clip2, op1);
    if ( !LineSameSide(clip1, clip2, op2, refp) ) {
	tmp = LineDist(clip1, clip2, op2);
	if (tmp > extra)
	    extra = tmp;
    }
    if (extra > 0) {
	clip1 = BP_ADD(refp, BP_SCALE(ut, min+extra));
	clip2 = BP_ADD(clip1, UT_90CW(ut));
    }
    intersects = IntersectLines(np1, &op1, &cop1, &clip1, &clip2);
    assert(intersects);
    intersects = IntersectLines(np2, &op2, &cop2, &clip1, &clip2);
    assert(intersects);
}

typedef struct joinparams {
    StrokeContext *c;
    SplineSet *cur;
    BasePoint sxy, oxy, ut_endlast;
    NibOffset  *noi;
    int is_ccw, was_ccw, is_right, bend_is_ccw;
} JoinParams;

static void BevelJoin(JoinParams *jpp) {
    SplineSetLineTo(jpp->cur, jpp->oxy);
}

static void NibJoin(JoinParams *jpp) {
    NibOffset now;
    SplinePoint *sp;

    CalcNibOffset(jpp->c, jpp->ut_endlast, jpp->is_right, &now, -1);
    sp = AddNibPortion(jpp->c->nibcorners, jpp->cur->last, &now, jpp->was_ccw,
                       jpp->noi, jpp->is_ccw, !jpp->bend_is_ccw);
    SplineStrokeAppendFixup(sp, jpp->sxy, jpp->noi, jpp->is_ccw);
    jpp->cur->last = sp;
}

static void DoubleBackJoin(JoinParams *jpp) {
    BasePoint refp, p1, p2;
    bigreal fsw, jlim, d1 = 0, d2 = 0;

    assert( jpp->c->join==lj_miterclip || jpp->c->join==lj_arcs );
    assert( RealWithin( BP_DOT(jpp->ut_endlast, jpp->noi->utanvec),
                        -1, COS_MARGIN) );
    refp = BP_ADD(jpp->sxy, jpp->c->pseudo_origin);
    fsw = FalseStrokeWidth(jpp->c, NormVec(BP_ADD(jpp->ut_endlast,
                                                  jpp->noi->utanvec)));
    jlim = CalcJoinLimit(jpp->c, fsw);
    CalcExtend(refp, jpp->ut_endlast, jpp->cur->last->me, jpp->oxy, jlim/2,
               &p1, &p2);
    SplineSetLineTo(jpp->cur, p1);
    SplineSetLineTo(jpp->cur, p2);
    SplineSetLineTo(jpp->cur, jpp->oxy);
}

static void RoundJoin(JoinParams *jpp) {
    BasePoint c, cut, ut1, ut2;
    bigreal alpha, B, C, E, mu, nu, B2AC, maj, min, tmp, tmp2;
    int intersects;

    c = BP_AVG(jpp->cur->last->me, jpp->oxy);
    cut = NormVec(BP_ADD(BP_REV(jpp->cur->last->me), jpp->oxy));
    alpha = BP_DIST(jpp->cur->last->me, jpp->oxy)/2;
    ut1 = BP_ROT(BP_REV(jpp->ut_endlast), UT_NEG(cut));
    ut2 = BP_ROT(jpp->noi->utanvec, UT_NEG(cut));
    BasePoint pp1 = BP_ROT(jpp->cur->last->me, UT_NEG(cut)), pp2 = BP_ROT(jpp->oxy, UT_NEG(cut));
    printf("p1: %lf,%lf p2: %lf,%lf pp1: %lf,%lf: pp2: %lf,%lf, c: %lf,%lf cut: %lf, alpha: %lf\n", jpp->cur->last->me.x, jpp->cur->last->me.y, jpp->oxy.x, jpp->oxy.y, pp1.x, pp1.y, pp2.x, pp2.y, c.x, c.y, atan2(cut.y, cut.x) * 180 / PI, alpha);
    mu = -ut1.x/ut1.y;
    nu = -ut2.x/ut2.y;
    B = 1 + pow(mu + nu, 2)/2;
    C = mu + nu;
    E = alpha * (mu - nu);
    B2AC = B*B-4*C;
    tmp = 2 * (E*E - B2AC*alpha*alpha);
    tmp2 = sqrt(pow(1-C,2) + B*B);
    printf("mu: %lf, nu:%lf, B: %lf, C: %lf, E: %lf, B2AC: %lf, tmp: %lf, tmp2: %lf\n", mu, nu, B, C, E, B2AC, tmp, tmp2);
    maj = sqrt(tmp * (1 + C + tmp2))/B2AC;
    min = sqrt(tmp * (1 + C - tmp2))/B2AC;
    printf("maj: %lf, min: %lf, angle: %lf\n", maj, min, atan2(cut.y, cut.x) * 180 / PI);
    BevelJoin(jpp);
} 

static void MiterJoin(JoinParams *jpp) {
    BasePoint ixy, refp, cow, coi, clip1, clip2, ut;
    int intersects;
    SplinePoint *sp;
    SplineSet *cur = jpp->cur;
    bigreal fsw, jlen, jlim;

    cow = BP_ADD(cur->last->me, jpp->ut_endlast);
    coi = BP_ADD(jpp->oxy, jpp->noi->utanvec);
    intersects = IntersectLines(&ixy, &cur->last->me, &cow, &coi, &jpp->oxy);
    assert(intersects); // Shouldn't be called with parallel tangents
    fsw = FalseStrokeWidth(jpp->c, NormVec(BP_ADD(jpp->ut_endlast,
                                                  jpp->noi->utanvec)));
    jlim = CalcJoinLimit(jpp->c, fsw);
    jlen = CalcJoinLength(fsw, jpp->ut_endlast, jpp->noi->utanvec);

    if ( jlen<=jlim ) {
	// Normal miter join
	SplineSetLineTo(cur, ixy);
	SplineSetLineTo(cur, jpp->oxy);
    } else {
	if ( jpp->c->join==lj_miter ) {
	    BevelJoin(jpp);
	    return;
	}
	refp = BP_ADD(jpp->sxy, jpp->c->pseudo_origin);
	ut = NormVec(BP_ADD(cur->last->me, BP_REV(jpp->oxy)));
	ut = jpp->bend_is_ccw ? UT_90CW(ut) : UT_90CCW(ut);
	clip1 = BP_ADD(refp, BP_SCALE(ut, jlim/2));
	clip2 = BP_ADD(clip1, UT_90CW(ut));
	if ( !LineSameSide(clip1, clip2, jpp->oxy, refp) ) {
	    // Don't trim past bevel
	    BevelJoin(jpp);
	    return;
	}
	if ( LineSameSide(clip1, clip2, ixy, refp) ) {
	    // Backup normal miter join (rounding problems)
	    SplineSetLineTo(cur, ixy);
	    SplineSetLineTo(cur, jpp->oxy);
	    return;
	}
	// Clipped miter join
	intersects = IntersectLines(&ixy, &cur->last->me, &cow, &clip1, &clip2);
	assert(intersects);
	SplineSetLineTo(cur, ixy);
	intersects = IntersectLines(&ixy, &jpp->oxy, &coi, &clip1, &clip2);
	assert(intersects);
	SplineSetLineTo(cur, ixy);
	SplineSetLineTo(cur, jpp->oxy);
    }
}

static int _HandleJoin(JoinParams *jpp) {
    SplineSet *cur = jpp->cur;
    BasePoint oxy = jpp->oxy;
    StrokeContext *c = jpp->c;
    int is_flat = false;
    bigreal costheta = BP_DOT(jpp->ut_endlast, jpp->noi->utanvec);

    if ( cur->first==NULL ) {
	// Create initial spline point
	cur->first = SplinePointCreate(oxy.x, oxy.y);
	cur->last = cur->first;
    } else if ( BPWITHIN(cur->last->me, oxy, INTERSPLINE_MARGIN) ) {
	// Close enough to just move the point
	SplineStrokeAppendFixup(cur->last, jpp->sxy, jpp->noi, jpp->is_ccw);
	is_flat = true;
    } else if ( RealWithin(costheta, 1, COS_MARGIN) ) {
        HandleFlat(jpp->cur, jpp->sxy, jpp->noi, jpp->is_ccw);
	is_flat = true;
    } else if ( !jpp->bend_is_ccw == !jpp->is_right ) {
	// Join is under nib, just connect for later removal
	BevelJoin(jpp);
    } else if ( c->join==lj_bevel ) {
	BevelJoin(jpp);
    } else if ( c->join==lj_round ) {
	if ( RealWithin(costheta, -1, COS_MARGIN) )
	    BevelJoin(jpp);
	else
	    RoundJoin(jpp);
    } else if ( c->join==lj_miter || c->join==lj_miterclip ) {
	if ( RealWithin(costheta, -1, COS_MARGIN) ) {
	    if ( c->join==lj_miter )
		BevelJoin(jpp);
	    else
		DoubleBackJoin(jpp);
	} else
	    MiterJoin(jpp);
    } else {
	if ( c->join!=lj_nib && c->join!=lj_inherited )
	    LogError( _("Warning: Unrecognized or unsupported join type, defaulting to 'nib'.\n") );
	NibJoin(jpp);
    } 
    return is_flat;
}

static int HandleJoin(StrokeContext *c, SplineSet *cur,
                      BasePoint sxy, NibOffset *noi, int is_ccw,
                      BasePoint ut_endlast, int was_ccw,
                      int is_right) {
    JoinParams jp = { c, cur, sxy, BP_UNINIT, ut_endlast, noi,
                      is_ccw, was_ccw, is_right, false };
    jp.oxy = BP_ADD(sxy, noi->off[is_ccw]);
    jp.bend_is_ccw = !JointBendsCW(ut_endlast, noi->utanvec);
    return _HandleJoin(&jp);
}

   
static void HandleCap(StrokeContext *c, SplineSet *cur, BasePoint sxy,
                      BasePoint ut, int is_ccw, int is_right) {
    NibOffset nof, not;
    BasePoint refp = BP_ADD(sxy, c->pseudo_origin), p1, p2, oxy;
    SplinePoint *sp;
    bigreal cw, cel;

    CalcNibOffset(c, ut, false, &nof, -1);
    CalcNibOffset(c, ut, true, &not, -1);
    oxy = BP_ADD(sxy, not.off[is_ccw]);

    if ( c->cap==lc_bevel ) {
	SplineSetLineTo(cur, oxy);
	return;
    }
    cw = BP_DIST(nof.off[is_ccw], not.off[is_ccw]);
    if ( c->cap==lc_butt || c->cap==lc_round ) {
	cel = CalcCapExtend(c, cw);
	CalcExtend(refp, BP_REV_IF(!is_right, ut),
	           cur->last->me, oxy, cel, &p1, &p2);
	if ( BPNEAR(cur->last->me, p1) )
	    cur->last->me = p1;
	else
	    SplineSetLineTo(cur, p1);
	if ( c->cap==lc_round ) {
	    SSAppendSemiCircle(cur, BP_DIST(p1, p2)/2, ut, !is_right);
	    assert( BPWITHIN(cur->last->me, p2, INTERSPLINE_MARGIN) );
	    cur->last->me = p2;
	} else {
	    SplineSetLineTo(cur, p2);
	}
	if ( !BPNEAR(oxy, p2) )
	    SplineSetLineTo(cur, oxy);
    } else {
	if ( c->cap!=lc_nib && c->cap!=lc_inherited )
	    LogError( _("Warning: Unrecognized or unsupported cap type, defaulting to 'nib'.\n") );
	sp = AddNibPortion(c->nibcorners, cur->last, &nof, is_ccw, &not,
	                   is_ccw, !is_right);
	// SplineStrokeAppendFixup(sp, sxy, &not, is_ccw);
	cur->last = sp;
    }
}

/* The main loop of the stroking algorithm.
 */
static SplineSet *OffsetSplineSet(SplineSet *ss, StrokeContext *c) {
    NibOffset no, no_last;
    Spline *s, *first=NULL;
    SplineSet *left=NULL, *right=NULL, *cur;
    SplinePoint *sp;
    BasePoint ut_ini = BP_UNINIT, ut_start, ut_mid, ut_end, ut_endlast;
    BasePoint sxy;
    bigreal last_t, t;
    int is_right, linear, curved, on_cusp;
    int is_ccw_ini, is_ccw_start, is_ccw_mid, was_ccw;
    int closed = ss->first->prev!=NULL;

    if ( (c->contour_was_ccw ? !c->remove_inner : !c->remove_outer) || !closed )
	left = chunkalloc(sizeof(SplineSet));
    if ( (c->contour_was_ccw ? !c->remove_outer : !c->remove_inner) || !closed )
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

	    if ( !HandleJoin(c, cur, sxy, &no, is_ccw_start, ut_endlast,
	                     was_ccw, is_right) )
		// Handle cusp value change here
		;
	    on_cusp = OffsetOnCuspAt(c, s, 0.0, &no, is_right, is_ccw_start);

	    printf("is_right=%d nci=%d is_ccw:%d nc.x=%lf nc.y=%lf, ut_start=%lf,%lf\n", is_right, no.nci[is_ccw_start], is_ccw_start, c->nibcorners[no.nci[is_ccw_start]].on_nib->me.x, c->nibcorners[no.nci[is_ccw_start]].on_nib->me.y, ut_start.x, ut_start.y);

	    // The path for this spline
	    if ( linear ) {
		sxy = SPLINEPVAL(s, 1.0);
		SplineSetLineTo(cur, BP_ADD(sxy, no.off[is_ccw_start]));
	    } else {
		t = 0.0;
		ut_mid = ut_start;
		is_ccw_mid = is_ccw_start;
		while ( t < 1.0 ) {
		    last_t = t;
		    t = SplineStrokeNextT(c, s, t, is_ccw_mid, &ut_mid, &curved,
		                          is_right, no.nci[is_ccw_mid]);
		    assert( t > last_t && t <= 1.0 );

		    printf("nci: %d, ut_mid=%lf,%lf, t=%.15lf, curved=%d\n", no.nci[is_ccw_mid], ut_mid.x, ut_mid.y, t, curved);

		    if ( curved )
			sp = TraceAndFitSpline(c, s, last_t, t, cur->last,
			                       no.nci[is_ccw_mid], is_right,
			                       on_cusp);
		    else
			sp = AppendCubicSplinePortion(s, last_t, t, cur->last);

		    cur->last = sp;

		    sxy = SPLINEPVAL(s, t);
		    CalcNibOffset(c, ut_mid, is_right, &no, no.nci[is_ccw_mid]);
		    is_ccw_mid = SplineTurningCCWAt(s, t);
		    SplineStrokeAppendFixup(cur->last, sxy, &no,
		                            on_cusp ? is_ccw_mid : !is_ccw_mid);

		    HandleFlat(cur, sxy, &no, is_ccw_mid);
		    on_cusp = OffsetOnCuspAt(c, s, t, &no, is_right, is_ccw_mid);
		}
	    }
	}
	ut_endlast = ut_end;
        was_ccw = SplineTurningCCWAt(s, 1.0);
    }
    if (    (left!=NULL && left->first==NULL) 
         || (right!=NULL && right->first==NULL) ) {
	// The path (presumably) had only zero-length splines
	LogError( _("Warning: No stroke output for contour\n") );
	assert(    (left==NULL || left->first==NULL)
	        && (right==NULL || right->first==NULL) );
	chunkfree(left, sizeof(SplineSet));
	chunkfree(right, sizeof(SplineSet));
	return NULL;
    }
    cur = NULL;
    if ( !closed ) {
	HandleCap(c, left, ss->last->me, ut_endlast, was_ccw, true);
	SplineSetReverse(right);
	left->next = right;
	right = NULL;
	SplineSetJoin(left, true, FIXUP_MARGIN, &closed);
	if ( !closed )
	     LogError( _("Warning: Contour end did not close\n") );
	else {
	    HandleCap(c, left, ss->first->me, ut_ini, is_ccw_ini, false);
	    SplineSetJoin(left, true, FIXUP_MARGIN, &closed);
	    if ( !closed )
		LogError( _("Warning: Contour start did not close\n") );
	    else if ( c->rmov==srmov_contour )
		left = SplineSetRemoveOverlap(NULL,left,over_remove);
	}
	cur = left;
	left = NULL;
    } else {
	// This can fail if the source contour is closed in a strange way
	if ( left!=NULL ) {
	    CalcNibOffset(c, ut_ini, false, &no, -1);
	    HandleJoin(c, left, ss->first->me, &no, is_ccw_ini, ut_endlast,
	               was_ccw, false);
            left = SplineSetJoin(left, true, INTERSPLINE_MARGIN, &closed);
	    if ( !closed )
		LogError( _("Warning: Left contour did not close\n") );
	    else if ( c->rmov==srmov_contour )
		left = SplineSetRemoveOverlap(NULL,left,over_remove);
	    cur = left;
	    left = NULL;
	}
	if ( right!=NULL ) {
	    CalcNibOffset(c, ut_ini, true, &no, -1);
	    HandleJoin(c, right, ss->first->me, &no, is_ccw_ini, ut_endlast,
	               was_ccw, true);
            right = SplineSetJoin(right, true, INTERSPLINE_MARGIN, &closed);
	    if ( !closed )
		LogError( _("Warning: Right contour did not close\n") );
	    else {
		SplineSetReverse(right);
		if ( c->rmov!=srmov_none )
		    // Need to do this for either srmov_contour or srmov_layer
		    right = SplineContourOuterCCWRemoveOverlap(right);
	    }
	    if ( cur != NULL ) {
		cur->next = right;
	    } else
		cur = right;
	    right = NULL;
	}
    }
    return cur;
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
    { { -1, 0 }, { -1, -CIRCOFF }, { -1, CIRCOFF } },
    { { 0 , 1 }, { -CIRCOFF, 1 }, { CIRCOFF, 1 } },
    { { 1, 0 }, { 1, CIRCOFF }, { 1, -CIRCOFF } },
    { { 0, -1 }, { CIRCOFF, -1 }, { -CIRCOFF, -1 } },
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

    c->contour_was_ccw = (   base->first->prev!=NULL
                          && !SplinePointListIsClockwise(base) );

    if ( base->first->next==NULL )
	ret = SinglePointStroke(base->first, c);
    else {
	if ( base->first->next->order2 ) {
	    base = SSPSApprox(ss);
	    copied = 1;
	}
	if ( c->contour_was_ccw )
	    SplineSetReverse(base);
	ret = OffsetSplineSet(base, c);
    }
    if ( order2 )
	ret = SplineSetsConvertOrder(ret,order2);
    if ( copied )
	SplinePointListFree(base);
    else if ( c->contour_was_ccw )
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
    bigreal sn = 0.0, co = 1.0, mr;
    DBounds b;
    real trans[6];
    struct simplifyinfo smpl = { sf_forcelines&sf_mergelines&sf_smoothcurves&sf_ignoreslopes, 0.25, 0.1, .005, 0, 0, 0 };

    if ( si->stroke_type==si_centerline )
	IError("centerline not handled");

    memset(&c,0,sizeof(c));
    c.join = si->join;
    c.cap  = si->cap;
    c.joinlimit = si->joinlimit;
    c.extendcap = si->extendcap;
    c.acctarget = si->accuracy_target;
    c.remove_inner = si->removeinternal;
    c.remove_outer = si->removeexternal;
    c.leave_users_center = si->leave_users_center;
    c.extrema = si->extrema;
    c.simplify = si->simplify;
    c.ecrelative = si->ecrelative;
    c.jlrelative = si->jlrelative;
    c.rmov = si->rmov;
    if ( si->minorradius!=0 )
	mr = si->minorradius;
    else
	mr = si->radius;
    //c.scaled_or_rotated = si->factor!=NULL;
    
    if ( c.acctarget<MIN_ACCURACY ) {
         LogError( _("Warning: Accuracy target %lf less than minimum %lf, using %lf instead\n"), c.acctarget, MIN_ACCURACY );
	 c.acctarget = MIN_ACCURACY;
    }
    smpl.err = c.acctarget;

    if ( si->penangle!=0 ) {
	sn = sin(si->penangle);
	co = cos(si->penangle);
    }
    trans[0] = co;
    trans[1] = sn;
    trans[2] = -sn;
    trans[3] = co;
    trans[4] = trans[5] = 0;

    if ( si->stroke_type==si_round || si->stroke_type==si_caligraphic ) {
	c.pentype = si->stroke_type==si_round ? pt_circle : pt_square;
	max_pc = 4;
	nibs = UnitShape(si->stroke_type==si_round ? 0 : -4);
	trans[0] *= si->radius;
	trans[1] *= si->radius;
	trans[2] *= mr;
	trans[3] *= mr;
    } else {
	c.pentype = pt_convex;
	max_pc = 20; // a guess
	nibs = SplinePointListCopy(si->nib);
	SplineSetFindBounds(nibs,&b);
	if ( !c.leave_users_center ) {
	    trans[4] = -(b.minx+b.maxx)/2;
	    trans[5] = -(b.miny+b.maxy)/2;
	} else {
	    c.pseudo_origin.x = b.minx+b.maxx/2;
	    c.pseudo_origin.y = b.miny+b.maxy/2;
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
	BuildNibCorners(&c.nibcorners, nib, &max_pc, &c.n);
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
    ff_progress_end_indicator();
}
