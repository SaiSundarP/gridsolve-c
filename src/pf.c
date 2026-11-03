#include "gridsolve.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static size_t bus_index(const GsCase *grid, int id)
{
    size_t i;
    for (i = 0; i < grid->nbus; ++i) if (grid->buses[i].id == id) return i;
    return SIZE_MAX;
}

static void injections(size_t n, const double complex *y, const double *vm,
                       const double *va, double *p, double *q)
{
    size_t i, j;
    for (i = 0; i < n; ++i) {
        p[i] = 0.0; q[i] = 0.0;
        for (j = 0; j < n; ++j) {
            double angle = va[i] - va[j];
            double g = creal(y[i*n+j]), b = cimag(y[i*n+j]);
            p[i] += vm[i] * vm[j] * (g*cos(angle) + b*sin(angle));
            q[i] += vm[i] * vm[j] * (g*sin(angle) - b*cos(angle));
        }
    }
}

GsStatus gs_solve(const GsCase *grid, const GsOptions *options, GsResult *out)
{
    GsOptions defaults = gs_default_options();
    const GsOptions *opt = options == NULL ? &defaults : options;
    double complex *y = NULL, *jac = NULL, *rhs = NULL;
    double *vm = NULL, *va = NULL, *p = NULL, *q = NULL, *qset = NULL, *history = NULL;
    int *types = NULL, *angle_col = NULL, *voltage_col = NULL;
    size_t n, maxvar, i, j, iter;
    GsStatus status;
    double mismatch = INFINITY;
    if (out == NULL) return GS_ERR_INVALID;
    *out = (GsResult){0};
    status = gs_case_validate(grid);
    if (status != GS_OK) return status;
    if (!isfinite(opt->tolerance) || opt->tolerance <= 0.0 || opt->max_iterations <= 0)
        return GS_ERR_INVALID;
    n = grid->nbus;
    if (n > SIZE_MAX/2 || 2*n > SIZE_MAX/(2*n)) return GS_ERR_NOMEM;
    if ((size_t)opt->max_iterations > SIZE_MAX/sizeof(*history)) return GS_ERR_NOMEM;
    maxvar = 2*n;
    status = gs_ybus(grid, &y);
    if (status != GS_OK) return status;
    vm = malloc(n*sizeof(*vm)); va = malloc(n*sizeof(*va));
    p = malloc(n*sizeof(*p)); q = malloc(n*sizeof(*q)); qset = malloc(n*sizeof(*qset));
    types = malloc(n*sizeof(*types)); angle_col = malloc(n*sizeof(*angle_col));
    voltage_col = malloc(n*sizeof(*voltage_col));
    jac = calloc(maxvar*maxvar, sizeof(*jac)); rhs = calloc(maxvar, sizeof(*rhs));
    history = malloc((size_t)opt->max_iterations*sizeof(*history));
    if (vm == NULL || va == NULL || p == NULL || q == NULL || qset == NULL || types == NULL ||
        angle_col == NULL || voltage_col == NULL || jac == NULL || rhs == NULL || history == NULL) {
        status = GS_ERR_NOMEM; goto cleanup;
    }
    for (i = 0; i < n; ++i) {
        types[i] = grid->buses[i].type;
        vm[i] = types[i] == GS_PQ ? 1.0 : grid->buses[i].vm;
        va[i] = opt->flat_start && types[i] != GS_SLACK ? 0.0 :
            grid->buses[i].va_deg * acos(-1.0) / 180.0;
        qset[i] = (grid->buses[i].qg - grid->buses[i].qd) / grid->base_mva;
    }
    for (iter = 0; iter < (size_t)opt->max_iterations; ++iter) {
        size_t na = 0, npq = 0, nvar;
        int switched = 0;
        for (i = 0; i < n; ++i) {
            angle_col[i] = -1; voltage_col[i] = -1;
            if (types[i] != GS_SLACK) angle_col[i] = (int)na++;
        }
        for (i = 0; i < n; ++i) if (types[i] == GS_PQ) voltage_col[i] = (int)(na + npq++);
        nvar = na + npq;
        injections(n, y, vm, va, p, q);
        if (opt->enforce_q_limits) {
            for (i = 0; i < n; ++i) if (types[i] == GS_PV) {
                double qgen = q[i] * grid->base_mva + grid->buses[i].qd;
                if (qgen < grid->buses[i].qmin || qgen > grid->buses[i].qmax) {
                    double limit = qgen < grid->buses[i].qmin ? grid->buses[i].qmin : grid->buses[i].qmax;
                    types[i] = GS_PQ;
                    vm[i] = grid->buses[i].vm;
                    qset[i] = (limit - grid->buses[i].qd) / grid->base_mva;
                    switched = 1;
                }
            }
        }
        if (switched) continue;
        mismatch = 0.0;
        memset(rhs, 0, maxvar*sizeof(*rhs));
        for (i = 0; i < n; ++i) {
            double p_spec = (grid->buses[i].pg - grid->buses[i].pd) / grid->base_mva;
            double q_spec = qset[i];
            double dp = p_spec - p[i];
            double dq = q_spec - q[i];
            if (types[i] != GS_SLACK) {
                rhs[(size_t)angle_col[i]] = dp;
                if (fabs(dp) > mismatch) mismatch = fabs(dp);
            }
            if (types[i] == GS_PQ) {
                rhs[(size_t)voltage_col[i]] = dq;
                if (fabs(dq) > mismatch) mismatch = fabs(dq);
            }
        }
        history[out->mismatch_count++] = mismatch;
        gs_log(GS_DEBUG, opt->log_level, "NR iteration %zu mismatch %.12e pu", iter, mismatch);
        if (mismatch <= opt->tolerance) {
            out->iterations = (int)iter;
            break;
        }
        if (nvar == 0) { status = GS_ERR_SINGULAR; goto cleanup; }
        memset(jac, 0, maxvar*maxvar*sizeof(*jac));
        {
            size_t p_row = 0, q_row = na;
            for (i = 0; i < n; ++i) if (types[i] != GS_SLACK) {
                for (j = 0; j < n; ++j) {
                    double angle = va[i] - va[j];
                    double g = creal(y[i*n+j]), b = cimag(y[i*n+j]);
                    if (i == j) {
                        jac[p_row*nvar+(size_t)angle_col[j]] = -q[i] - cimag(y[i*n+i])*vm[i]*vm[i];
                        if (types[i] == GS_PQ)
                            jac[p_row*nvar+(size_t)voltage_col[j]] = p[i]/vm[i] + creal(y[i*n+i])*vm[i];
                        if (types[i] == GS_PQ) {
                            jac[q_row*nvar+(size_t)angle_col[j]] = p[i] - creal(y[i*n+i])*vm[i]*vm[i];
                            jac[q_row*nvar+(size_t)voltage_col[j]] = q[i]/vm[i] - cimag(y[i*n+i])*vm[i];
                        }
                    } else {
                        if (angle_col[j] >= 0)
                            jac[p_row*nvar+(size_t)angle_col[j]] = vm[i]*vm[j]*(g*sin(angle)-b*cos(angle));
                        if (types[j] == GS_PQ)
                            jac[p_row*nvar+(size_t)voltage_col[j]] = vm[i]*(g*cos(angle)+b*sin(angle));
                        if (types[i] == GS_PQ && angle_col[j] >= 0)
                            jac[q_row*nvar+(size_t)angle_col[j]] = -vm[i]*vm[j]*(g*cos(angle)+b*sin(angle));
                        if (types[i] == GS_PQ && types[j] == GS_PQ)
                            jac[q_row*nvar+(size_t)voltage_col[j]] = vm[i]*(g*sin(angle)-b*cos(angle));
                    }
                }
                ++p_row;
                if (types[i] == GS_PQ) ++q_row;
            }
        }
        status = gs_lu_solve(nvar, jac, rhs);
        if (status != GS_OK) goto cleanup;
        for (i = 0; i < n; ++i) {
            if (angle_col[i] >= 0) va[i] += creal(rhs[(size_t)angle_col[i]]);
            if (voltage_col[i] >= 0) vm[i] += creal(rhs[(size_t)voltage_col[i]]);
            if (!isfinite(vm[i]) || vm[i] <= 0.0) { status = GS_ERR_NOCONV; goto cleanup; }
        }
    }
    if (iter == (size_t)opt->max_iterations && mismatch > opt->tolerance) {
        status = GS_ERR_NOCONV; goto cleanup;
    }
    injections(n, y, vm, va, p, q);
    out->vm = malloc(n*sizeof(*out->vm)); out->va_rad = malloc(n*sizeof(*out->va_rad));
    out->p_inj = malloc(n*sizeof(*out->p_inj)); out->q_inj = malloc(n*sizeof(*out->q_inj));
    if (out->vm == NULL || out->va_rad == NULL || out->p_inj == NULL || out->q_inj == NULL) {
        gs_result_free(out); status = GS_ERR_NOMEM; goto cleanup;
    }
    memcpy(out->vm, vm, n*sizeof(*vm)); memcpy(out->va_rad, va, n*sizeof(*va));
    memcpy(out->p_inj, p, n*sizeof(*p)); memcpy(out->q_inj, q, n*sizeof(*q));
    out->mismatch_history = history; history = NULL;
    out->nbus = n; out->nbranch = grid->nbranch; out->mismatch = mismatch;
    {
        GsFlow *flows = calloc(grid->nbranch == 0 ? 1 : grid->nbranch, sizeof(*flows));
        if (flows == NULL) { gs_result_free(out); status = GS_ERR_NOMEM; goto cleanup; }
        status = gs_branch_flows(grid, out, flows);
        if (status == GS_OK) for (i = 0; i < grid->nbranch; ++i) if (grid->branches[i].status) {
            out->total_loss_p += flows[i].p_from + flows[i].p_to;
            out->total_loss_q += flows[i].q_from + flows[i].q_to;
        }
        free(flows);
        if (status != GS_OK) { gs_result_free(out); goto cleanup; }
    }
    status = GS_OK;
cleanup:
    free(y); free(vm); free(va); free(p); free(q); free(qset); free(types);
    free(angle_col); free(voltage_col); free(jac); free(rhs); free(history);
    if (status != GS_OK) gs_result_free(out);
    return status;
}

GsStatus gs_branch_flows(const GsCase *grid, const GsResult *result, GsFlow *flows)
{
    size_t i;
    if (grid == NULL || result == NULL || flows == NULL || result->vm == NULL ||
        result->va_rad == NULL || result->nbus != grid->nbus || result->nbranch != grid->nbranch)
        return GS_ERR_INVALID;
    for (i = 0; i < grid->nbranch; ++i) {
        const GsBranch *br = &grid->branches[i];
        size_t f = bus_index(grid, br->from), t = bus_index(grid, br->to);
        double shift;
        double complex tap, ys, yff, yft, ytf, ytt, vf, vt, iff, itt;
        if (f == SIZE_MAX || t == SIZE_MAX) return GS_ERR_INVALID;
        flows[i] = (GsFlow){0};
        if (!br->status) continue;
        shift = br->shift_deg*acos(-1.0)/180.0;
        tap = br->tap*(cos(shift)+I*sin(shift));
        ys = 1.0/(br->r+I*br->x);
        ytt = ys+I*br->b/2.0; yff = ytt/(tap*conj(tap));
        yft = -ys/conj(tap); ytf = -ys/tap;
        vf = result->vm[f]*(cos(result->va_rad[f])+I*sin(result->va_rad[f]));
        vt = result->vm[t]*(cos(result->va_rad[t])+I*sin(result->va_rad[t]));
        iff = yff*vf+yft*vt; itt = ytf*vf+ytt*vt;
        flows[i].p_from = creal(vf*conj(iff))*grid->base_mva;
        flows[i].q_from = cimag(vf*conj(iff))*grid->base_mva;
        flows[i].p_to = creal(vt*conj(itt))*grid->base_mva;
        flows[i].q_to = cimag(vt*conj(itt))*grid->base_mva;
    }
    return GS_OK;
}