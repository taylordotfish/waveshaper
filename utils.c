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

#include "utils.h"

#define NUM_LOWPASS_STAGES 3

static int compare_points(const void *p1_ptr, const void *p2_ptr) {
    const Point *point1 = (const Point *)p1_ptr;
    const Point *point2 = (const Point *)p2_ptr;

    int cmp = double_cmp(point1->x, point2->x);
    if (cmp != 0) return cmp;

    if (point1->index < point2->index) return 1;
    if (point1->index > point2->index) return -1;
    return 0;
}

void sort_points(Point *points, size_t count) {
    qsort(points, count, sizeof(points[0]), &compare_points);
}

size_t make_points_unique(Point *points, size_t count) {
    size_t offset = 0;
    for (size_t i = 1; i < count;) {
        points[i] = points[i + offset];
        if (double_cmp(points[i].x, points[i - 1].x) == 0) {
            offset++;
            count--;
        } else {
            i++;
        }
    }
    return count;
}

size_t sort_points_unique(Point *points, size_t count) {
    sort_points(points, count);
    return make_points_unique(points, count);
}

void init_lowpass(LowpassParams *params) {
    params->gt = init_iir_stage(
        IIR_STAGE_LOWPASS, NUM_LOWPASS_STAGES, 3, 2
    );
    params->iirf = init_iirf_t(params->gt);
}

void update_lowpass_coefficients(LowpassParams *params, float fc) {
    chebyshev(
        params->iirf, params->gt, NUM_LOWPASS_STAGES * 2, IIR_STAGE_LOWPASS,
        fc, 0.5f
    );
}

void run_lowpass(LowpassParams *params, const float *input, float *output,
                 uint32_t n_samples) {
    iir_process_buffer_ns_5(
        params->iirf, params->gt, input, output, n_samples
    );
}

void reset_lowpass(LowpassParams *params) {
    reset_iirf_t(params->iirf, params->gt, params->gt->availst);
}

void free_lowpass(LowpassParams *params) {
    free_iirf_t(params->iirf, params->gt);
    free_iir_stage(params->gt);
    params->gt = NULL;
    params->iirf = NULL;
}

void upsample_raw(const float *input, uint32_t n_samples,
                  uint_fast8_t multiplier, float *output) {
    for (size_t pos = 0, base = 0; pos < n_samples; pos++) {
        float val = input[pos];
        for (size_t offset = 0; offset < multiplier; offset++) {
            output[base + offset] = val;
        }
        base += multiplier;
    }
}

void upsample(const float *input, uint32_t n_samples,
              uint_fast8_t multiplier, float *output, LowpassParams *lp) {
    float *upsamples = malloc(n_samples * multiplier * sizeof(float));
    upsample_raw(input, n_samples, multiplier, upsamples);
    run_lowpass(lp, upsamples, output, n_samples * multiplier);
    free(upsamples);
}

void downsample_raw(const float *input, uint32_t n_samples,
                    uint_fast8_t multiplier, float *output) {
    for (size_t pos = 0, new_pos = 0; pos < n_samples; pos += multiplier) {
        output[new_pos] = input[pos];
        new_pos++;
    }
}

void downsample(const float *input, uint32_t n_samples,
                uint_fast8_t multiplier, float *output, LowpassParams *lp) {
    float *lpsamples = malloc(n_samples * sizeof(float));
    run_lowpass(lp, input, lpsamples, n_samples);
    downsample_raw(lpsamples, n_samples, multiplier, output);
    free(lpsamples);
}
