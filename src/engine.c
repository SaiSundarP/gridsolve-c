#include "gridsolve.h"

#include <math.h>
#include <stdlib.h>

GsStatus gs_connectivity(const GsCase *grid)
{
    return gs_case_validate(grid);
}

static int rank_results(const void *left, const void *right)
{
    const GsContingencyResult *a = left, *b = right;
    if (a->status != GS_OK && b->status == GS_OK) return -1;
    if (a->status == GS_OK && b->status != GS_OK) return 1;
    if (a->severity < b->severity) return 1;
    if (a->severity > b->severity) return -1;
    return a->branch_index - b->branch_index;
}

static GsStatus assess(const GsCase *grid, const GsResult *solution,
                       GsContingencyResult *item)
{
    GsFlow *flows;
    size_t i;
    GsStatus status;
    flows = calloc(grid->nbranch == 0 ? 1 : grid->nbranch, sizeof(*flows));
    if (flows == NULL) return GS_ERR_NOMEM;
    status = gs_branch_flows(grid, solution, flows);
    if (status != GS_OK) { free(flows); return status; }
    for (i = 0; i < grid->nbus; ++i) {
        double low = grid->buses[i].vmin-solution->vm[i];
        double high = solution->vm[i]-grid->buses[i].vmax;
        if (low > 0.0) { ++item->voltage_violations; item->severity += low; }
        if (high > 0.0) { ++item->voltage_violations; item->severity += high; }
    }
    for (i = 0; i < grid->nbranch; ++i) if (grid->branches[i].status && grid->branches[i].rate_a > 0.0) {
        double from = hypot(flows[i].p_from, flows[i].q_from);
        double to = hypot(flows[i].p_to, flows[i].q_to);
        double maximum = fmax(from, to);
        if (maximum > grid->branches[i].rate_a) {
            ++item->overloads;
            item->severity += maximum/grid->branches[i].rate_a-1.0;
        }
    }
    free(flows);
    return GS_OK;
}

GsStatus gs_n1(const GsCase *grid, const GsOptions *options,
               GsContingencyResult **out, size_t *count)
{
    GsCase outage = {0};
    GsOptions flat_options = options == NULL ? gs_default_options() : *options;
    GsContingencyResult *results = NULL;
    size_t i, n = 0, k = 0;
    GsStatus status;
    if (out == NULL || count == NULL) return GS_ERR_INVALID;
    *out = NULL; *count = 0;
    status = gs_case_validate(grid);
    if (status != GS_OK) return status;
    flat_options.flat_start = 1;
    for (i = 0; i < grid->nbranch; ++i) if (grid->branches[i].status) ++n;
    results = calloc(n == 0 ? 1 : n, sizeof(*results));
    if (results == NULL) return GS_ERR_NOMEM;
    status = gs_case_clone(grid, &outage);
    if (status != GS_OK) { free(results); return status; }
    for (i = 0; i < grid->nbranch; ++i) if (grid->branches[i].status) {
        GsResult solution = {0};
        GsContingencyResult *item = &results[k++];
        item->branch_index = (int)i;
        outage.branches[i].status = 0;
        status = gs_connectivity(&outage);
        if (status == GS_ERR_ISLAND) item->status = GS_ERR_ISLAND;
        else if (status != GS_OK) { gs_case_free(&outage); free(results); return status; }
        else {
            item->status = gs_solve(&outage, &flat_options, &solution);
            if (item->status == GS_OK) {
                status = assess(&outage, &solution, item);
                if (status != GS_OK) {
                    gs_result_free(&solution); gs_case_free(&outage); free(results); return status;
                }
            }
        }
        gs_result_free(&solution);
        outage.branches[i].status = 1;
    }
    gs_case_free(&outage);
    qsort(results, n, sizeof(*results), rank_results);
    *out = results; *count = n;
    return GS_OK;
}

GsStatus gs_sweep(const GsCase *grid, const GsOptions *options,
                  double start, double stop, double step,
                  double *first_violation, double *collapse_point)
{
    GsCase scaled = {0};
    GsOptions flat_options = options == NULL ? gs_default_options() : *options;
    double factor;
    GsStatus status;
    if (first_violation == NULL || collapse_point == NULL || !isfinite(start) ||
        !isfinite(stop) || !isfinite(step) || start <= 0.0 || stop < start || step <= 0.0 ||
        start+step <= start || (stop-start)/step > 1000000.0)
        return GS_ERR_INVALID;
    *first_violation = NAN; *collapse_point = NAN;
    status = gs_case_validate(grid);
    if (status != GS_OK) return status;
    flat_options.flat_start = 1;
    status = gs_case_clone(grid, &scaled);
    if (status != GS_OK) return status;
    for (factor = start; factor <= stop + step*1e-9; factor += step) {
        GsResult solution = {0};
        GsContingencyResult assessment = {0};
        size_t i;
        for (i = 0; i < scaled.nbus; ++i) {
            scaled.buses[i].pd = grid->buses[i].pd*factor;
            scaled.buses[i].qd = grid->buses[i].qd*factor;
        }
        status = gs_solve(&scaled, &flat_options, &solution);
        if (status != GS_OK) {
            *collapse_point = factor;
            gs_result_free(&solution);
            break;
        }
        status = assess(&scaled, &solution, &assessment);
        gs_result_free(&solution);
        if (status != GS_OK) break;
        if (assessment.voltage_violations > 0 || assessment.overloads > 0) {
            if (isnan(*first_violation)) *first_violation = factor;
        }
    }
    gs_case_free(&scaled);
    return status == GS_ERR_NOCONV || status == GS_ERR_SINGULAR ? GS_OK : status;
}