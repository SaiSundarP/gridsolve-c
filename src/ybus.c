#include "gridsolve.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>

static size_t index_of(const GsCase *grid, int id)
{
    size_t i;
    for (i = 0; i < grid->nbus; ++i) if (grid->buses[i].id == id) return i;
    return SIZE_MAX;
}

GsStatus gs_ybus(const GsCase *grid, double complex **out)
{
    size_t i, n;
    double complex *y;
    if (out == NULL) return GS_ERR_INVALID;
    *out = NULL;
    if (grid == NULL || grid->nbus == 0 || grid->buses == NULL)
        return GS_ERR_INVALID;
    n = grid->nbus;
    {
        GsStatus status = gs_case_validate(grid);
        if (status != GS_OK) return status;
    }
    if (n > SIZE_MAX / n || n * n > SIZE_MAX / sizeof(*y)) return GS_ERR_NOMEM;
    y = calloc(n * n, sizeof(*y));
    if (y == NULL) return GS_ERR_NOMEM;
    for (i = 0; i < n; ++i) y[i*n+i] += grid->buses[i].gs + I * grid->buses[i].bs;
    for (i = 0; i < grid->nbranch; ++i) {
        const GsBranch *br = &grid->branches[i];
        size_t f, t;
        double complex ys, tap;
        double shift;
        if (!br->status) continue;
        f = index_of(grid, br->from); t = index_of(grid, br->to);
        if (f == SIZE_MAX || t == SIZE_MAX) { free(y); return GS_ERR_INVALID; }
        ys = 1.0 / (br->r + I * br->x);
        shift = br->shift_deg * acos(-1.0) / 180.0;
        tap = br->tap * (cos(shift) + I * sin(shift));
        y[f*n+f] += (ys + I * br->b/2.0) / (tap * conj(tap));
        y[f*n+t] += -ys / conj(tap);
        y[t*n+f] += -ys / tap;
        y[t*n+t] += ys + I * br->b/2.0;
    }
    *out = y;
    return GS_OK;
}