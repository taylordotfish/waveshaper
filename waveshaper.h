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

#ifndef WAVESHAPER_H
#define WAVESHAPER_H

#include "cspline.h"
#include "ports.h"
#include "utils.h"
#include <lv2/lv2plug.in/ns/lv2core/lv2.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

// Constants
enum {
    NUM_POINTS = NUM_USER_POINTS + 3,
    NUM_SPLINES = NUM_POINTS - 1,
    PORTS_PER_POINT = 6,
};

typedef struct {
    const float *xs[NUM_POINTS];
    const float *ys[NUM_POINTS];
    const float *types[NUM_SPLINES];
    const float *tensions[NUM_SPLINES];
    const float *tangent1_vals[NUM_SPLINES];
    const float *tangent2_vals[NUM_SPLINES];
    const float *display_plot;
    const float *oversample;
    const float *input_gain;
    const float *output_gain;
    const float *monitor_levels;
    float *input_level;
    float *input_clip;
    float *output_level;
    float *output_0db;
    const float *asymmetric;

    float old_xs[NUM_POINTS];
    float old_ys[NUM_POINTS];
    float old_types[NUM_SPLINES];
    float old_tensions[NUM_SPLINES];
    float old_tangent1_vals[NUM_SPLINES];
    float old_tangent2_vals[NUM_SPLINES];
    bool old_display_plot;
    uint_fast8_t old_oversample_multiplier;
    bool old_monitor_levels;
    bool old_asymmetric;

    uint32_t sample_rate;
    const float *inputs[2];
    float *outputs[2];

    SplineParams splines[NUM_POINTS + 1];
    size_t num_splines;
    FILE *gnuplot;

    LowpassParams lp1[2];
    LowpassParams lp2[2];

    float *oversample_buf;
    size_t oversample_buf_size;

    uint32_t reset_input_clip_after;
    uint32_t reset_output_0db_after;
} Waveshaper;

int gnuplot_printf(Waveshaper *ws, const char *format, ...);
void start_gnuplot(Waveshaper *ws);
void update_gnuplot(Waveshaper *ws);
void stop_gnuplot(Waveshaper *ws);
void close_gnuplot(Waveshaper *ws);
bool needs_update(Waveshaper *ws);
void update_splines(Waveshaper *ws);
void update_caches(Waveshaper *ws);
double get_curve_value(Waveshaper *ws, double x);
void update_lp_filters(Waveshaper *ws, uint_fast8_t multiplier);
void reset_lp_filters(Waveshaper *ws);

LV2_SYMBOL_EXPORT
const LV2_Descriptor *lv2_descriptor(uint32_t index);

#endif
