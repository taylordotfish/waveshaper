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

#define _POSIX_C_SOURCE 200112L

#include "waveshaper.h"
#include <math.h>
#include <signal.h>
#include <stdarg.h>
#include <string.h>
#include <unistd.h>

#pragma GCC diagnostic ignored "-Wunused-parameter"

#define NAME "Waveshaper (Fish)"
#define URI "https://taylor.fish/plugins/waveshaper"
#define CLIP_TIMEOUT 2  // seconds
#define MONITOR_MIN_DB -40

static float f_zero = 0.0f;
static float f_half = 0.5f;
static float f_one = 1.0f;

int gnuplot_printf(Waveshaper *ws, const char *format, ...) {
    va_list args;
    va_start(args, format);
    int status = vfprintf(ws->gnuplot, format, args);
    va_end(args);
    if (status < 0) close_gnuplot(ws);
    return status;
}

void start_gnuplot(Waveshaper *ws) {
    if (ws->gnuplot != NULL) return;
    FILE *gnuplot = popen("LD_LIBRARY_PATH= gnuplot", "w");
    setbuf(gnuplot, NULL);
    ws->gnuplot = gnuplot;
    if (gnuplot_printf(ws, "unset key\n") < 0) return;
    if (gnuplot_printf(ws, "h00(t) = 2*t**3 - 3*t**2 + 1\n") < 0) return;
    if (gnuplot_printf(ws, "h10(t) = t**3 - 2*t**2 + t\n") < 0) return;
    if (gnuplot_printf(ws, "h01(t) = -2*t**3 + 3*t**2\n") < 0) return;
    if (gnuplot_printf(ws, "h11(t) = t**3 - t**2\n") < 0) return;

    if (gnuplot_printf(
        ws, "t(x, x1, x2) = (x - x1) / (x2 - x1)\n"
        "a(v, x1, x2, y1, y2) = (y2 - y1) / (x2 - x1) * v\n"
    ) < 0) return;

    if (gnuplot_printf(
        ws, "pc(x, x1, x2, y1, y2, m1, m2) = "
        "h00(t(x, x1, x2)) * y1 + h10(t(x, x1, x2)) * "
        "(x2 - x1) * a(m1 + 1, x1, x2, y1, y2) + "
        "h01(t(x, x1, x2)) * y2 + h11(t(x, x1, x2)) * "
        "(x2 - x1) * a(m2 + 1, x1, x2, y1, y2)\n"
    ) < 0) return;

    if (gnuplot_printf(
        ws, "pq(x, x1, x2, y1, y2, T) = "
        "(x2 - x1) * (-a(T, x1, x2, y1, y2) * t(x, x1, x2)**2 + "
        "a(T, x1, x2, y1, y2) * t(x, x1, x2)) - "
        "y1 * t(x, x1, x2) + y2 * t(x, x1, x2) + y1\n"
    ) < 0) return;

    if (gnuplot_printf(
        ws, "pl(x, x1, x2, y1, y2) = "
        "(y2 - y1) / (x2 - x1) * (x - x1) + y1\n"
    ) < 0) return;

    update_gnuplot(ws);
}

void update_gnuplot(Waveshaper *ws) {
    if (ws->gnuplot == NULL) return;
    if (gnuplot_printf(ws, "plot [0:1] [0:1] sample ") < 0) return;
    for (size_t i = 0; i < ws->num_splines; i++) {
        double x1 = ws->splines[i].x1;
        double x2 = ws->splines[i].x2;
        double y1 = ws->splines[i].y1;
        double y2 = ws->splines[i].y2;
        double tension = ws->splines[i].tension;
        double tangent1 = ws->splines[i].tangent1;
        double tangent2 = ws->splines[i].tangent2;

        if (gnuplot_printf(ws, "[%f:%f] ", x1, x2) < 0) return;
        int status = 0;
        switch (ws->splines[i].type) {
            case SPLINE_CUBIC:
                status = gnuplot_printf(
                    ws, "pc(x, %f, %f, %f, %f, %f, %f), ", x1, x2, y1, y2,
                    tangent1, tangent2
                );
                break;

            case SPLINE_QUADRATIC:
                status = gnuplot_printf(
                    ws, "pq(x, %f, %f, %f, %f, %f), ", x1, x2, y1, y2, tension
                );
                break;

            default:
                status = gnuplot_printf(
                    ws, "pl(x, %f, %f, %f, %f), ", x1, x2, y1, y2
                );
                break;
        }
        if (status < 0) return;
    }

    if (gnuplot_printf(ws, "\"-\" with points lt 7\n") < 0) return;
    for (size_t i = 0; i < ws->num_splines; i++) {
        if (i == 0) {
            double x1 = ws->splines[i].x1;
            double y1 = ws->splines[i].y1;
            if (gnuplot_printf(ws, "%f %f\n", x1, y1) < 0) return;
        }

        double x2 = ws->splines[i].x2;
        double y2 = ws->splines[i].y2;
        if (gnuplot_printf(ws, "%f %f\n", x2, y2) < 0) return;
    }
    if (gnuplot_printf(ws, "e\n") < 0) return;
}

void stop_gnuplot(Waveshaper *ws) {
    if (ws->gnuplot == NULL) return;
    if (gnuplot_printf(ws, "\nquit\n") < 0) return;
    close_gnuplot(ws);
}

void close_gnuplot(Waveshaper *ws) {
    pclose(ws->gnuplot);
    ws->gnuplot = NULL;
}

bool needs_update(Waveshaper *ws) {
    if ((*ws->asymmetric > 0) != ws->old_asymmetric) {
        return true;
    }

    for (size_t i = 0; i < NUM_POINTS; i++) {
        if (float_cmp(*ws->xs[i], ws->old_xs[i]) != 0) {
            return true;
        }
        if (float_cmp(*ws->ys[i], ws->old_ys[i]) != 0) {
            return true;
        }
    }

    for (size_t i = 0; i < NUM_SPLINES; i++) {
        if (float_cmp(*ws->types[i], ws->old_types[i]) != 0) {
            return true;
        }

        SplineType type = (SplineType)*ws->types[i];
        if (type < 0 || type >= NUM_SPLINETYPES) type = SPLINE_LINEAR;

        if (type == SPLINE_QUADRATIC) {
            if (float_cmp(*ws->tensions[i], ws->old_tensions[i]) != 0) {
                return true;
            }
        }

        if (type == SPLINE_CUBIC) {
            double tangent1 = *ws->tangent1_vals[i];
            double tangent2 = *ws->tangent2_vals[i];
            if (float_cmp(tangent1, ws->old_tangent1_vals[i]) != 0) {
                return true;
            }
            if (float_cmp(tangent2, ws->old_tangent2_vals[i]) != 0) {
                return true;
            }
        }
    }

    return false;
}

void update_splines(Waveshaper *ws) {
    bool asymmetric = *ws->asymmetric > 0;
    if (asymmetric) {
        ws->xs[NUM_USER_POINTS + 1] = &f_half;
        ws->ys[NUM_USER_POINTS + 1] = &f_half;
    } else {
        ws->xs[NUM_USER_POINTS + 1] = &f_zero;
        ws->ys[NUM_USER_POINTS + 1] = &f_zero;
    }

    Point points[NUM_POINTS];
    for (size_t i = 0; i < NUM_POINTS; i++) {
        points[i].x = *ws->xs[i];
        points[i].y = *ws->ys[i];
        points[i].index = i;
    }

    size_t num_points = sort_points_unique(points, NUM_POINTS);
    for (size_t i = 0; i < num_points - 1; i++) {
        size_t index = points[i].index;
        SplineType type = (SplineType)*ws->types[index];
        if (type < 0 || type >= NUM_SPLINETYPES) type = SPLINE_LINEAR;

        ws->splines[i].x1 = points[i].x;
        ws->splines[i].x2 = points[i + 1].x;
        ws->splines[i].y1 = points[i].y;
        ws->splines[i].y2 = points[i + 1].y;
        ws->splines[i].tension = *ws->tensions[index];
        ws->splines[i].tangent1 = *ws->tangent1_vals[index];
        ws->splines[i].tangent2 = *ws->tangent2_vals[index];
        ws->splines[i].type = type;
    }
    ws->num_splines = num_points - 1;
}

void update_caches(Waveshaper *ws) {
    ws->old_asymmetric = *ws->asymmetric > 0;

    for (size_t i = 0; i < NUM_POINTS; i++) {
        ws->old_xs[i] = *ws->xs[i];
        ws->old_ys[i] = *ws->ys[i];
    }

    for (size_t i = 0; i < NUM_SPLINES; i++) {
        ws->old_types[i] = *ws->types[i];
        ws->old_tensions[i] = *ws->tensions[i];
        ws->old_tangent1_vals[i] = *ws->tangent1_vals[i];
        ws->old_tangent2_vals[i] = *ws->tangent2_vals[i];
    }
}

double get_curve_value(Waveshaper *ws, double x) {
    for (size_t i = 0; i < ws->num_splines; i++) {
        double x1 = ws->splines[i].x1;
        double x2 = ws->splines[i].x2;
        bool spline_matches = (
            (i == 0 && x <= x2) ||
            (i == ws->num_splines - 1) ||
            (x <= x2 && x >= x1)
        );
        if (spline_matches) {
            return spline_from_params(x, &ws->splines[i]);
        }
    }
    return 0;
}

void update_lp_filters(Waveshaper *ws, uint_fast8_t multiplier) {
    float fc = 1.0f / (2 * multiplier);
    for (size_t i = 0; i < 2; i++) {
        update_lowpass_coefficients(&ws->lp1[i], fc);
        update_lowpass_coefficients(&ws->lp2[i], fc);
    }
}

void reset_lp_filters(Waveshaper *ws) {
    for (size_t i = 0; i < 2; i++) {
        reset_lowpass(&ws->lp1[i]);
        reset_lowpass(&ws->lp2[i]);
    }
}

void free_lp_filters(Waveshaper *ws) {
    for (size_t i = 0; i < 2; i++) {
        free_lowpass(&ws->lp1[i]);
        free_lowpass(&ws->lp2[i]);
    }
}

static LV2_Handle instantiate(
        const LV2_Descriptor *descriptor, double rate, const char *bundle_path,
        const LV2_Feature * const *features) {
    signal(SIGPIPE, SIG_IGN);
    Waveshaper *ws = calloc(1, sizeof(Waveshaper));

    ws->sample_rate = rate;
    ws->oversample_buf_size = rate / 2;
    ws->oversample_buf = malloc(ws->oversample_buf_size * sizeof(float));
    for (size_t i = 0; i < 2; i++) {
        init_lowpass(&ws->lp1[i]);
        init_lowpass(&ws->lp2[i]);
    }

    ws->xs[NUM_USER_POINTS] = &f_zero;  // Asym start point
    ws->xs[NUM_USER_POINTS + 1] = &f_zero;  // Origin
    ws->ys[NUM_USER_POINTS + 1] = &f_zero;  // Origin
    ws->xs[NUM_USER_POINTS + 2] = &f_one;  // End point

    // This causes monitors to be reset during the first call to run().
    ws->old_monitor_levels = true;
    return (LV2_Handle)ws;
}

static void connect_port(LV2_Handle instance, uint32_t port, void *data) {
    Waveshaper *ws = (Waveshaper *)instance;

    bool is_point_setting = true;
    size_t index = 0;
    uint32_t port_type = 0;

    if (port >= MIN_POINT_PORT && port <= MAX_POINT_PORT) {
        index = (port - MIN_POINT_PORT) / PORTS_PER_POINT;
        port_type = (port - MIN_POINT_PORT) % PORTS_PER_POINT;
    } else if (port >= MIN_ASYM_PORT && port <= MAX_ASYM_PORT) {
        index = NUM_USER_POINTS;
        port_type = (port - MIN_ASYM_PORT) + 1;  // Starts with Y
    } else if (port >= MIN_ORIGIN_PORT && port <= MAX_ORIGIN_PORT) {
        index = NUM_USER_POINTS + 1;
        port_type = (port - MIN_ORIGIN_PORT) + 2;  // Starts with tension
    } else if (port == PORT_END_POINT_Y) {
        index = NUM_USER_POINTS + 2;
        port_type = 1;  // Y
    } else {
        is_point_setting = false;
    }

    if (is_point_setting) {
        switch (port_type) {
            case 0:
                ws->xs[index] = data;
                break;
            case 1:
                ws->ys[index] = data;
                break;
            case 2:
                ws->tensions[index] = data;
                break;
            case 3:
                ws->tangent1_vals[index] = data;
                break;
            case 4:
                ws->tangent2_vals[index] = data;
                break;
            case 5:
                ws->types[index] = data;
                break;
        }
        return;
    }

    switch (port) {
        case PORT_INPUT_L:
            ws->inputs[0] = data;
            break;
        case PORT_OUTPUT_L:
            ws->outputs[0] = data;
            break;
        case PORT_INPUT_R:
            ws->inputs[1] = data;
            break;
        case PORT_OUTPUT_R:
            ws->outputs[1] = data;
            break;
        case PORT_DISPLAY_PLOT:
            ws->display_plot = data;
            break;
        case PORT_OVERSAMPLE:
            ws->oversample = data;
            break;
        case PORT_INPUT_GAIN:
            ws->input_gain = data;
            break;
        case PORT_OUTPUT_GAIN:
            ws->output_gain = data;
            break;
        case PORT_MONITOR_LEVELS:
            ws->monitor_levels = data;
            break;
        case PORT_INPUT_LEVEL:
            ws->input_level = data;
            break;
        case PORT_INPUT_CLIP:
            ws->input_clip = data;
            break;
        case PORT_OUTPUT_LEVEL:
            ws->output_level = data;
            break;
        case PORT_OUTPUT_0DB:
            ws->output_0db = data;
            break;
        case PORT_ASYMMETRIC:
            ws->asymmetric = data;
            break;
    }
}

static void activate(LV2_Handle instance) {
    Waveshaper *ws = (Waveshaper *)instance;

    memset(ws->old_xs, 0, sizeof(ws->old_xs));
    memset(ws->old_ys, 0, sizeof(ws->old_ys));
    memset(ws->old_types, 0, sizeof(ws->old_types));
    memset(ws->old_tensions, 0, sizeof(ws->old_tensions));
    memset(ws->old_tangent1_vals, 0, sizeof(ws->old_tangent1_vals));
    memset(ws->old_tangent2_vals, 0, sizeof(ws->old_tangent2_vals));
    ws->old_display_plot = false;
    ws->reset_input_clip_after = 0;
    ws->reset_output_0db_after = 0;

    memset(ws->splines, 0, sizeof(ws->splines));
    ws->splines[0].x1 = 0;
    ws->splines[0].y1 = 0;
    ws->splines[0].x2 = 1;
    ws->splines[0].y2 = 1;
    ws->splines[0].type = SPLINE_LINEAR;
    ws->num_splines = 1;

    stop_gnuplot(ws);
    reset_lp_filters(ws);

    if (ws->input_level != NULL) *ws->input_level = 0;
    if (ws->input_clip != NULL) *ws->input_clip = 0;
    if (ws->output_level != NULL) *ws->output_level = 0;
    if (ws->output_0db != NULL) *ws->output_0db = 0;
}

static void run(LV2_Handle instance, uint32_t n_samples_total) {
    if (n_samples_total == 0) return;
    Waveshaper *ws = (Waveshaper *)instance;
    const float * const *inputs = ws->inputs;
    float * const *outputs = ws->outputs;
    bool display_plot = *ws->display_plot > 0;
    bool monitor = *ws->monitor_levels > 0;
    bool asymmetric = *ws->asymmetric > 0;

    uint_fast8_t multiplier = 0;
    if (*ws->oversample > 1 && *ws->oversample <= UINT8_MAX) {
        multiplier = (uint_fast8_t)*ws->oversample;
        if (multiplier != ws->old_oversample_multiplier) {
            update_lp_filters(ws, multiplier);
            ws->old_oversample_multiplier = multiplier;
        }
    }

    if (needs_update(ws)) {
        update_splines(ws);
        update_caches(ws);
        update_gnuplot(ws);
    }

    if (display_plot != ws->old_display_plot) {
        ws->old_display_plot = display_plot;
        if (display_plot) {
            start_gnuplot(ws);
        } else {
            stop_gnuplot(ws);
        }
    }

    uint32_t n_samples_max = n_samples_total;
    if (multiplier > 0) {
        if (n_samples_max * multiplier > ws->oversample_buf_size) {
            n_samples_max = ws->oversample_buf_size / multiplier;
        }
    }

    float max_in_val = 0;
    float max_out_val = 0;

    for (size_t port = 0; port < 2; port++) {
        const float *input = inputs[port];
        float *output = outputs[port];

        uint32_t n_samples = n_samples_max;
        uint32_t n_processed = 0;

        while (n_processed < n_samples_total) {
            if (n_samples > n_samples_total - n_processed) {
                n_samples = n_samples_total - n_processed;
            }

            uint32_t n_upsamples = n_samples;
            const float *upsamples = input;
            float *upsamples_out = output;

            if (multiplier > 0) {
                n_upsamples *= multiplier;
                upsample(
                    input, n_samples, multiplier, ws->oversample_buf,
                    &ws->lp1[port]
                );
                upsamples = ws->oversample_buf;
                upsamples_out = ws->oversample_buf;
            }

            for (uint32_t pos = 0; pos < n_upsamples; pos++) {
                float in_val = upsamples[pos];
                if (in_val == 0) {
                    upsamples_out[pos] = 0;
                    continue;
                }

                in_val *= db_to_amplitude(*ws->input_gain);
                int in_sign = in_val < 0 ? -1 : 1;

                if (monitor && fabs(in_val) > max_in_val) {
                    max_in_val = fabs(in_val);
                }

                if (in_val >= 1 || in_val <= -1) in_val = in_sign;
                if (asymmetric) {
                    in_val = (in_val + 1) / 2;
                } else {
                    in_val *= in_sign;
                }

                double out_val = get_curve_value(ws, in_val);
                if (asymmetric) {
                    out_val = (out_val * 2) - 1;
                } else {
                    out_val *= in_sign;
                }

                if (out_val >= 1) {
                    out_val = 1;
                } else if (out_val <= -1) {
                    out_val = -1;
                }

                upsamples_out[pos] = out_val;
            }

            if (multiplier > 0) {
                downsample(
                    upsamples_out, n_upsamples, multiplier, output,
                    &ws->lp2[port]
                );
            }

            for (uint32_t pos = 0; pos < n_samples; pos++) {
                float out_val = output[pos];
                out_val *= db_to_amplitude(*ws->output_gain);
                if (monitor) {
                    if (out_val > max_out_val || -out_val > max_out_val) {
                        max_out_val = fabs(out_val);
                    }
                }
                output[pos] = out_val;
            }

            input += n_samples;
            output += n_samples;
            n_processed += n_samples;
        }
    }

    if (monitor != ws->old_monitor_levels) {
        if (!monitor) {
            if (ws->input_clip != NULL) *ws->input_clip = 0;
            if (ws->output_0db != NULL) *ws->output_0db = 0;
            if (ws->input_level != NULL) *ws->input_level = MONITOR_MIN_DB;
            if (ws->output_level != NULL) *ws->output_level = MONITOR_MIN_DB;
            ws->reset_input_clip_after = 0;
            ws->reset_output_0db_after = 0;
        }
        ws->old_monitor_levels = monitor;
    }

    if (monitor) {
        if (n_samples_total >= ws->reset_input_clip_after) {
            ws->reset_input_clip_after = 0;
            if (ws->input_clip != NULL) *ws->input_clip = 0;
        } else if (ws->reset_input_clip_after > 0) {
            ws->reset_input_clip_after -= n_samples_total;
        }

        if (n_samples_total >= ws->reset_output_0db_after) {
            ws->reset_output_0db_after = 0;
            if (ws->output_0db != NULL) *ws->output_0db = 0;
        } else if (ws->reset_output_0db_after > 0) {
            ws->reset_output_0db_after -= n_samples_total;
        }

        if (max_in_val > 1) {
            if (ws->input_level != NULL) *ws->input_level = 0;
            if (ws->input_clip != NULL) *ws->input_clip = 1;
            ws->reset_input_clip_after = CLIP_TIMEOUT * ws->sample_rate;
        } else if (ws->input_level != NULL) {
            *ws->input_level = amplitude_to_db(max_in_val, MONITOR_MIN_DB);
        }

        if (max_out_val >= 1) {
            if (ws->output_level != NULL) *ws->output_level = 0;
            if (ws->output_0db != NULL) *ws->output_0db = 1;
            ws->reset_output_0db_after = CLIP_TIMEOUT * ws->sample_rate;
        } else if (ws->output_level != NULL) {
            *ws->output_level = amplitude_to_db(max_out_val, MONITOR_MIN_DB);
        }
    }
}

static void deactivate(LV2_Handle instance) {
    Waveshaper *ws = (Waveshaper *)instance;
    stop_gnuplot(ws);
}

static void cleanup(LV2_Handle instance) {
    Waveshaper *ws = (Waveshaper *)instance;
    free_lp_filters(ws);
    free(ws->oversample_buf);
    free(ws);
}

static const void *extension_data(const char *uri) {
    return NULL;
}

static const LV2_Descriptor descriptor = {
    URI,
    instantiate,
    connect_port,
    activate,
    run,
    deactivate,
    cleanup,
    extension_data
};

LV2_SYMBOL_EXPORT
const LV2_Descriptor *lv2_descriptor(uint32_t index) {
    switch(index) {
        case 0:
            return &descriptor;
        default:
            return NULL;
    }
}
