/*
 * Copyright (C) 2017 taylor.fish <contact@taylor.fish>
 *
 * This file is part of Fish Waveshaper.
 *
 * Fish Waveshaper is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Fish Waveshaper is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Fish Waveshaper.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CSPLINE_H
#define CSPLINE_H

typedef enum {
    SPLINE_CUBIC,
    SPLINE_QUADRATIC,
    SPLINE_LINEAR,
    NUM_SPLINETYPES
} SplineType;

typedef struct {
    double x1;
    double x2;
    double y1;
    double y2;
    double tension;
    double tangent1;
    double tangent2;
    SplineType type;
} SplineParams;

double cspline(double x, double x1, double x2, double y1, double y2,
               double tangent1, double tangent2);

double cspline_quad(double x, double x1, double x2, double y1, double y2,
                    double tension);

double linear_spline(double x, double x1, double x2, double y1, double y2);
double spline_from_params(double x, SplineParams *params);

#endif
