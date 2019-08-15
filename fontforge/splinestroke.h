#ifndef FONTFORGE_SPLINESTROKE_H
#define FONTFORGE_SPLINESTROKE_H

#include "splinefont.h"

enum ShapeType {
	Shape_Convex,
	Shape_Concave,
	Shape_PointOnEdge,
	Shape_TooFewPoints,
	Shape_Line
};

extern enum ShapeType NibIsValid(SplineSet *);
extern SplinePoint *AppendCubicSplinePortion(Spline *s, bigreal t_start,
                                             bigreal t_end, SplinePoint *start);
extern SplinePoint *AppendCubicSplineSetPortion(Spline *s, bigreal t_start,
                                                Spline *s_end, bigreal t_end,
                                                SplinePoint *dst_start,
                                                int backward);
extern SplineSet *SplineSetStroke(SplineSet *ss, StrokeInfo *si, int order2);
extern SplineSet *UnitShape(int n);
extern void FVStrokeItScript(void *_fv, StrokeInfo *si, int pointless_argument);

#endif /* FONTFORGE_SPLINESTROKE_H */
