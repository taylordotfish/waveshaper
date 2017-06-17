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

#include "cspline.h"

static inline double get_t(double x, double x1, double x2) {
    return (x - x1) / (x2 - x1);
}

double cspline(double x, double x1, double x2, double y1, double y2,
               double tangent1, double tangent2) {
    double slope = (y2 - y1) / (x2 - x1);
    double abs_tangent1 = slope * (tangent1 + 1);
    double abs_tangent2 = slope * (tangent2 + 1);
    double t = get_t(x, x1, x2);
    double t_2 = t * t;
    double t_3 = t_2 * t;
    return (
        (2 * t_3 - 3 * t_2 + 1) * y1 +
        (t_3 - 2 * t_2 + t) * (x2 - x1) * abs_tangent1 +
        (-2 * t_3 + 3 * t_2) * y2 +
        (t_3 - t_2) * (x2 - x1) * abs_tangent2
    );
}

double cspline_quad(double x, double x1, double x2, double y1, double y2,
                    double tension) {
    double slope = (y2 - y1) / (x2 - x1);
    double abs_tension = slope * tension;
    double t = get_t(x, x1, x2);
    double t_2 = t * t;
    return (
        (x2 - x1) * (-abs_tension * t_2 + abs_tension * t) -
        y1 * t + y2 * t + y1
    );
}

double linear_spline(double x, double x1, double x2, double y1, double y2) {
    double slope = (y2 - y1) / (x2 - x1);
    return slope * (x - x1) + y1;
}

double spline_from_params(double x, SplineParams *params) {
    double x1 = params->x1;
    double x2 = params->x2;
    double y1 = params->y1;
    double y2 = params->y2;

    switch (params->type) {
        case SPLINE_CUBIC:
            return cspline(
                x, x1, x2, y1, y2, params->tangent1, params->tangent2
            );
        case SPLINE_QUADRATIC:
            return cspline_quad(x, x1, x2, y1, y2, params->tension);
        default:
            return linear_spline(x, x1, x2, y1, y2);
    }
}
