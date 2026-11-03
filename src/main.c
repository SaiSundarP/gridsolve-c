#include "gridsolve.h"

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int parse_double(const char *text, double *value)
{
    char *end;
    errno = 0;
    *value = strtod(text, &end);
    return errno == 0 && end != text && *end == '\0' && isfinite(*value);
}

static int parse_int(const char *text, int *value)
{
    char *end;
    long parsed;
    errno = 0;
    parsed = strtol(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' || parsed < 1 || parsed > 1000000) return 0;
    *value = (int)parsed;
    return 1;
}

static void usage(const char *program)
{
    (void)fprintf(stderr, "Usage: %s pf|n1|sweep <case> [--tol value] [--maxit n] [--qlim] [--csv] [--step value]\n", program);
}

static GsStatus print_pf(const GsCase *grid, const GsResult *result, int csv)
{
    size_t i;
    GsFlow *flows = calloc(grid->nbranch == 0 ? 1 : grid->nbranch, sizeof(*flows));
    GsStatus status;
    if (flows == NULL) return GS_ERR_NOMEM;
    status = gs_branch_flows(grid, result, flows);
    if (status != GS_OK) { free(flows); return status; }
    if (csv) {
        (void)puts("record,id_or_from,to,vm_pu,va_deg,pg_mw,qg_mvar,p_from_mw,q_from_mvar,p_to_mw,q_to_mvar");
        for (i = 0; i < grid->nbus; ++i)
            (void)printf("bus,%d,,%.10f,%.10f,%.8f,%.8f,,,,\n", grid->buses[i].id,
                         result->vm[i], result->va_rad[i]*180.0/acos(-1.0),
                         result->p_inj[i]*grid->base_mva+grid->buses[i].pd,
                         result->q_inj[i]*grid->base_mva+grid->buses[i].qd);
        for (i = 0; i < grid->nbranch; ++i)
            (void)printf("branch,%d,%d,,,,,%.8f,%.8f,%.8f,%.8f\n", grid->branches[i].from,
                         grid->branches[i].to, flows[i].p_from, flows[i].q_from,
                         flows[i].p_to, flows[i].q_to);
    } else {
        (void)printf("Converged in %d iterations; mismatch %.3e pu\n", result->iterations, result->mismatch);
        (void)printf("Bus    Vm(pu)    Va(deg)       Pg(MW)     Qg(MVAr)\n");
        for (i = 0; i < grid->nbus; ++i)
            (void)printf("%-6d %8.5f %10.5f %12.4f %12.4f\n", grid->buses[i].id,
                         result->vm[i], result->va_rad[i]*180.0/acos(-1.0),
                         result->p_inj[i]*grid->base_mva+grid->buses[i].pd,
                         result->q_inj[i]*grid->base_mva+grid->buses[i].qd);
        (void)printf("Total losses: %.6f MW, %.6f MVAr\n", result->total_loss_p, result->total_loss_q);
        (void)puts("From  To       Pfrom       Qfrom         Pto         Qto");
        for (i = 0; i < grid->nbranch; ++i)
            (void)printf("%-5d %-5d %11.4f %11.4f %11.4f %11.4f\n", grid->branches[i].from,
                         grid->branches[i].to, flows[i].p_from, flows[i].q_from,
                         flows[i].p_to, flows[i].q_to);
    }
    free(flows);
    return GS_OK;
}

int main(int argc, char **argv)
{
    const char *mode, *case_path;
    GsCase grid = {0};
    GsOptions options = gs_default_options();
    double step = 0.05;
    int csv = 0, i;
    GsStatus status;
    if (argc < 3) { usage(argv[0]); return EXIT_FAILURE; }
    mode = argv[1]; case_path = argv[2];
    for (i = 3; i < argc; ++i) {
        if (strcmp(argv[i], "--qlim") == 0) options.enforce_q_limits = 1;
        else if (strcmp(argv[i], "--csv") == 0) csv = 1;
        else if ((strcmp(argv[i], "--tol") == 0 || strcmp(argv[i], "--step") == 0 ||
                  strcmp(argv[i], "--maxit") == 0) && i+1 < argc) {
            const char *option = argv[i++];
            if (strcmp(option, "--maxit") == 0) {
                if (!parse_int(argv[i], &options.max_iterations)) { usage(argv[0]); return EXIT_FAILURE; }
            } else {
                double value;
                if (!parse_double(argv[i], &value) || value <= 0.0) { usage(argv[0]); return EXIT_FAILURE; }
                if (strcmp(option, "--tol") == 0) options.tolerance = value;
                else step = value;
            }
        } else { usage(argv[0]); return EXIT_FAILURE; }
    }
    status = gs_case_load(case_path, &grid);
    if (status != GS_OK) {
        gs_log(GS_ERROR, options.log_level, "%s: %s", case_path, gs_status_string(status));
        return (int)status;
    }
    if (strcmp(mode, "pf") == 0) {
        GsResult result = {0};
        status = gs_solve(&grid, &options, &result);
        if (status == GS_OK) status = print_pf(&grid, &result, csv);
        gs_result_free(&result);
    } else if (strcmp(mode, "n1") == 0) {
        GsContingencyResult *results = NULL;
        size_t count = 0, j;
        status = gs_n1(&grid, &options, &results, &count);
        if (status == GS_OK) {
            if (csv) (void)puts("rank,from,to,status,voltage_violations,overloads,severity");
            else (void)puts("Rank From To Status VoltageViolations Overloads Severity");
            for (j = 0; j < count; ++j) {
                const GsContingencyResult *item = &results[j];
                const GsBranch *br = &grid.branches[item->branch_index];
                if (csv) (void)printf("%zu,%d,%d,%s,%d,%d,%.8f\n", j+1, br->from, br->to,
                    gs_status_string(item->status), item->voltage_violations, item->overloads, item->severity);
                else (void)printf("%4zu %4d %2d %-34s %4d %8d %.8f\n", j+1, br->from, br->to,
                    gs_status_string(item->status), item->voltage_violations, item->overloads, item->severity);
            }
        }
        free(results);
    } else if (strcmp(mode, "sweep") == 0) {
        double first, collapse;
        status = gs_sweep(&grid, &options, 0.5, 1.5, step, &first, &collapse);
        if (status == GS_OK) {
            if (csv) {
                (void)puts("first_violation_scale,collapse_scale");
                (void)printf("%.6f,%.6f\n", first, collapse);
            } else {
                if (isnan(first)) (void)puts("First violation scale: not reached");
                else (void)printf("First violation scale: %.6f\n", first);
                if (isnan(collapse)) (void)puts("Collapse scale: not reached");
                else (void)printf("Collapse scale: %.6f\n", collapse);
            }
        }
    } else { usage(argv[0]); status = GS_ERR_INVALID; }
    gs_case_free(&grid);
    if (status != GS_OK) {
        gs_log(GS_ERROR, options.log_level, "%s", gs_status_string(status));
        return (int)status;
    }
    return EXIT_SUCCESS;
}