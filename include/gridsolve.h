#ifndef GRIDSOLVE_H
#define GRIDSOLVE_H

#include <complex.h>
#include <stddef.h>

typedef enum {
    GS_OK = 0,
    GS_ERR_IO = 1,
    GS_ERR_PARSE = 2,
    GS_ERR_INVALID = 3,
    GS_ERR_NOMEM = 4,
    GS_ERR_SINGULAR = 5,
    GS_ERR_NOCONV = 6,
    GS_ERR_ISLAND = 7
} GsStatus;

typedef enum { GS_PQ = 1, GS_PV = 2, GS_SLACK = 3 } GsBusType;
typedef enum { GS_ERROR = 0, GS_WARN = 1, GS_INFO = 2, GS_DEBUG = 3 } GsLogLevel;

typedef struct {
    int id, type;
    double vm, va_deg, pd, qd, pg, qg, gs, bs, vmin, vmax, qmin, qmax;
} GsBus;

typedef struct {
    int from, to, status;
    double r, x, b, tap, shift_deg, rate_a;
} GsBranch;

typedef struct {
    double base_mva;
    size_t nbus, nbranch;
    GsBus *buses;
    GsBranch *branches;
} GsCase;

typedef struct {
    double tolerance;
    int max_iterations;
    int enforce_q_limits;
    int flat_start;
    int log_level;
} GsOptions;

typedef struct {
    int iterations;
    double mismatch;
    size_t nbus, nbranch;
    double *vm, *va_rad, *p_inj, *q_inj;
    double *mismatch_history;
    size_t mismatch_count;
    double total_loss_p, total_loss_q;
} GsResult;

typedef struct {
    double p_from, q_from, p_to, q_to;
} GsFlow;

typedef struct {
    int branch_index;
    GsStatus status;
    int voltage_violations, overloads;
    double severity;
} GsContingencyResult;

void gs_log(int level, int configured_level, const char *format, ...);
const char *gs_status_string(GsStatus status);
GsOptions gs_default_options(void);
GsStatus gs_case_load(const char *path, GsCase *out);
void gs_case_free(GsCase *grid);
GsStatus gs_case_validate(const GsCase *grid);
GsStatus gs_case_clone(const GsCase *source, GsCase *out);
GsStatus gs_ybus(const GsCase *grid, double complex **out);
GsStatus gs_lu_solve(size_t n, double complex *matrix, double complex *rhs);
GsStatus gs_solve(const GsCase *grid, const GsOptions *options, GsResult *out);
void gs_result_free(GsResult *result);
GsStatus gs_branch_flows(const GsCase *grid, const GsResult *result, GsFlow *flows);
GsStatus gs_connectivity(const GsCase *grid);
GsStatus gs_n1(const GsCase *grid, const GsOptions *options,
               GsContingencyResult **out, size_t *count);
GsStatus gs_sweep(const GsCase *grid, const GsOptions *options,
                  double start, double stop, double step,
                  double *first_violation, double *collapse_point);

#endif