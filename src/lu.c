#include "gridsolve.h"

#include <math.h>
#include <stdint.h>

GsStatus gs_lu_solve(size_t n, double complex *matrix, double complex *rhs)
{
    size_t i, j, k;
    double scale = 0.0, min_pivot = INFINITY;
    if (n == 0 || matrix == NULL || rhs == NULL || n > SIZE_MAX / n) return GS_ERR_INVALID;
    for (i = 0; i < n*n; ++i) {
        double v = cabs(matrix[i]);
        if (!isfinite(v)) return GS_ERR_INVALID;
        if (v > scale) scale = v;
    }
    for (i = 0; i < n; ++i)
        if (!isfinite(creal(rhs[i])) || !isfinite(cimag(rhs[i]))) return GS_ERR_INVALID;
    if (scale == 0.0) return GS_ERR_SINGULAR;
    for (k = 0; k < n; ++k) {
        size_t pivot = k;
        double largest = cabs(matrix[k*n+k]);
        for (i = k + 1; i < n; ++i) {
            double candidate = cabs(matrix[i*n+k]);
            if (candidate > largest) { largest = candidate; pivot = i; }
        }
        if (!isfinite(largest) || largest <= scale * 1e-12) return GS_ERR_SINGULAR;
        if (largest < min_pivot) min_pivot = largest;
        if (pivot != k) {
            double complex tmp;
            for (j = k; j < n; ++j) { tmp = matrix[k*n+j]; matrix[k*n+j] = matrix[pivot*n+j]; matrix[pivot*n+j] = tmp; }
            tmp = rhs[k]; rhs[k] = rhs[pivot]; rhs[pivot] = tmp;
        }
        for (i = k + 1; i < n; ++i) {
            double complex factor = matrix[i*n+k] / matrix[k*n+k];
            matrix[i*n+k] = 0.0;
            for (j = k + 1; j < n; ++j) matrix[i*n+j] -= factor * matrix[k*n+j];
            rhs[i] -= factor * rhs[k];
        }
    }
    if (min_pivot / scale <= 1e-12) return GS_ERR_SINGULAR;
    for (i = n; i-- > 0;) {
        for (j = i + 1; j < n; ++j) rhs[i] -= matrix[i*n+j] * rhs[j];
        rhs[i] /= matrix[i*n+i];
        if (!isfinite(creal(rhs[i])) || !isfinite(cimag(rhs[i]))) return GS_ERR_SINGULAR;
    }
    return GS_OK;
}