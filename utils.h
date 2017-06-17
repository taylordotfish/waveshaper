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

#ifndef UTILS_H
#define UTILS_H

#include "iir.h"
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#define FLOAT_CMP_PRECISION 0.00001
#define DOUBLE_CMP_PRECISION 0.00001

typedef struct {
    double x;
    double y;
    size_t index;
} Point;

typedef struct {
    iir_stage_t *gt;
    iirf_t *iirf;
} LowpassParams;

static inline int float_cmp(float first, float second) {
    float diff = first - second;
    if (diff < -FLOAT_CMP_PRECISION) return -1;
    if (diff > FLOAT_CMP_PRECISION) return 1;
    return 0;
}

static inline int double_cmp(double first, double second) {
    double diff = first - second;
    if (diff < -DOUBLE_CMP_PRECISION) return -1;
    if (diff > DOUBLE_CMP_PRECISION) return 1;
    return 0;
}

static inline float db_multiplier(float db) {
    return powf(10.0f, db * 0.05f);
}

void sort_points(Point *points, size_t count);
size_t make_points_unique(Point *points, size_t count);
size_t sort_points_unique(Point *points, size_t count);

void init_lowpass(LowpassParams *params);
void update_lowpass_coefficients(LowpassParams *params, float fc);
void run_lowpass(LowpassParams *params, const float *input, float *output,
                 uint32_t n_samples);
void reset_lowpass(LowpassParams *params);
void free_lowpass(LowpassParams *params);

void upsample_raw(const float *input, uint32_t n_samples,
                  uint_fast8_t multiplier, float *output);

void upsample(const float *input, uint32_t n_samples,
              uint_fast8_t multiplier, float *output, LowpassParams *lp);

void downsample_raw(const float *input, uint32_t n_samples,
                    uint_fast8_t multiplier, float *output);

void downsample(const float *input, uint32_t n_samples,
                uint_fast8_t multiplier, float *output, LowpassParams *lp);

#endif
