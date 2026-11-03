#include "gridsolve.h"

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

void gs_case_free(GsCase *grid)
{
    if (grid == NULL) return;
    free(grid->buses);
    free(grid->branches);
    *grid = (GsCase){0};
}

static int bus_index(const GsCase *grid, int id)
{
    size_t i;
    for (i = 0; i < grid->nbus; ++i) if (grid->buses[i].id == id) return (int)i;
    return -1;
}

GsStatus gs_case_validate(const GsCase *grid)
{
    size_t i, j;
    int slack_count = 0;
    unsigned char *seen = NULL;
    size_t *queue = NULL;
    size_t head = 0, tail = 0;
    if (grid == NULL || grid->buses == NULL || grid->nbus == 0 ||
        !isfinite(grid->base_mva) || grid->base_mva <= 0.0 ||
        (grid->nbranch > 0 && grid->branches == NULL)) return GS_ERR_INVALID;
    for (i = 0; i < grid->nbus; ++i) {
        const GsBus *b = &grid->buses[i];
        const double values[] = {b->vm,b->va_deg,b->pd,b->qd,b->pg,b->qg,b->gs,b->bs,
                                 b->vmin,b->vmax,b->qmin,b->qmax};
        if (b->id <= 0 || (b->type != GS_PQ && b->type != GS_PV && b->type != GS_SLACK))
            return GS_ERR_INVALID;
        for (j = 0; j < sizeof(values)/sizeof(values[0]); ++j)
            if (!isfinite(values[j])) return GS_ERR_INVALID;
        if (b->vm <= 0.0 || b->vmin <= 0.0 || b->vmax < b->vmin ||
            b->qmin > b->qmax) return GS_ERR_INVALID;
        if (b->type == GS_SLACK) ++slack_count;
        for (j = 0; j < i; ++j) if (grid->buses[j].id == b->id) return GS_ERR_INVALID;
    }
    if (slack_count != 1) return GS_ERR_INVALID;
    for (i = 0; i < grid->nbranch; ++i) {
        const GsBranch *br = &grid->branches[i];
        if ((br->status != 0 && br->status != 1) || bus_index(grid, br->from) < 0 ||
            bus_index(grid, br->to) < 0 || br->from == br->to ||
            !isfinite(br->r) || !isfinite(br->x) || !isfinite(br->b) ||
            !isfinite(br->tap) || !isfinite(br->shift_deg) || !isfinite(br->rate_a) ||
            br->r < 0.0 || br->x < 0.0 || (br->r == 0.0 && br->x == 0.0) ||
            br->tap <= 0.0 || br->rate_a < 0.0) return GS_ERR_INVALID;
    }
    seen = calloc(grid->nbus, sizeof(*seen));
    queue = malloc(grid->nbus * sizeof(*queue));
    if (seen == NULL || queue == NULL) { free(seen); free(queue); return GS_ERR_NOMEM; }
    seen[0] = 1U;
    queue[tail++] = 0U;
    while (head < tail) {
        size_t here = queue[head++];
        for (i = 0; i < grid->nbranch; ++i) {
            const GsBranch *br = &grid->branches[i];
            int a, z;
            size_t next;
            if (!br->status) continue;
            a = bus_index(grid, br->from); z = bus_index(grid, br->to);
            if (a < 0 || z < 0) { free(seen); free(queue); return GS_ERR_INVALID; }
            if ((size_t)a == here) next = (size_t)z;
            else if ((size_t)z == here) next = (size_t)a;
            else continue;
            if (!seen[next]) { seen[next] = 1U; queue[tail++] = next; }
        }
    }
    for (i = 0; i < grid->nbus; ++i) if (!seen[i]) {
        free(seen); free(queue); return GS_ERR_ISLAND;
    }
    free(seen); free(queue);
    return GS_OK;
}

GsStatus gs_case_clone(const GsCase *source, GsCase *out)
{
    if (source == NULL || out == NULL) return GS_ERR_INVALID;
    if ((source->nbus > 0 && source->buses == NULL) ||
        (source->nbranch > 0 && source->branches == NULL) ||
        source->nbus > SIZE_MAX/sizeof(*source->buses) ||
        source->nbranch > SIZE_MAX/sizeof(*source->branches)) return GS_ERR_INVALID;
    *out = (GsCase){0};
    out->base_mva = source->base_mva;
    out->nbus = source->nbus; out->nbranch = source->nbranch;
    if (out->nbus > 0) {
        out->buses = malloc(out->nbus * sizeof(*out->buses));
        if (out->buses == NULL) goto nomem;
        memcpy(out->buses, source->buses, out->nbus * sizeof(*out->buses));
    }
    if (out->nbranch > 0) {
        out->branches = malloc(out->nbranch * sizeof(*out->branches));
        if (out->branches == NULL) goto nomem;
        memcpy(out->branches, source->branches, out->nbranch * sizeof(*out->branches));
    }
    return GS_OK;
nomem:
    gs_case_free(out);
    return GS_ERR_NOMEM;
}

GsStatus gs_case_load(const char *path, GsCase *out)
{
    FILE *file;
    char line[1024];
    size_t line_number = 0, bus_capacity = 0, branch_capacity = 0;
    GsCase grid = {0};
    int have_base = 0;
    GsStatus status = GS_OK;
    if (path == NULL || out == NULL) return GS_ERR_INVALID;
    *out = (GsCase){0};
    file = fopen(path, "r");
    if (file == NULL) return GS_ERR_IO;
    while (fgets(line, sizeof(line), file) != NULL) {
        char *p = line;
        ++line_number;
        while (isspace((unsigned char)*p)) ++p;
        if (*p == '\0' || *p == '#' || *p == '\n') continue;
        if (strncmp(p, "BASEMVA", 7) == 0) {
            char extra;
            if (sscanf(p, "BASEMVA %lf %c", &grid.base_mva, &extra) != 1 || have_base) {
                status = GS_ERR_PARSE; break;
            }
            have_base = 1;
        } else if (strncmp(p, "BUS", 3) == 0) {
            GsBus b;
            char extra;
            if (sscanf(p, "BUS %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %c",
                       &b.id,&b.type,&b.vm,&b.va_deg,&b.pd,&b.qd,&b.pg,&b.qg,&b.gs,&b.bs,
                       &b.vmin,&b.vmax,&b.qmin,&b.qmax,&extra) != 14) { status = GS_ERR_PARSE; break; }
            if (grid.nbus == bus_capacity) {
                size_t next = bus_capacity == 0 ? 8 : bus_capacity * 2;
                GsBus *grown;
                if (next < bus_capacity || next > SIZE_MAX / sizeof(*grown)) { status = GS_ERR_NOMEM; break; }
                grown = realloc(grid.buses, next * sizeof(*grown));
                if (grown == NULL) { status = GS_ERR_NOMEM; break; }
                grid.buses = grown; bus_capacity = next;
            }
            grid.buses[grid.nbus++] = b;
        } else if (strncmp(p, "BRANCH", 6) == 0) {
            GsBranch b;
            char extra;
            if (sscanf(p, "BRANCH %d %d %lf %lf %lf %lf %lf %lf %d %c",
                       &b.from,&b.to,&b.r,&b.x,&b.b,&b.tap,&b.shift_deg,&b.rate_a,&b.status,&extra) != 9) {
                status = GS_ERR_PARSE; break;
            }
            if (grid.nbranch == branch_capacity) {
                size_t next = branch_capacity == 0 ? 8 : branch_capacity * 2;
                GsBranch *grown;
                if (next < branch_capacity || next > SIZE_MAX / sizeof(*grown)) { status = GS_ERR_NOMEM; break; }
                grown = realloc(grid.branches, next * sizeof(*grown));
                if (grown == NULL) { status = GS_ERR_NOMEM; break; }
                grid.branches = grown; branch_capacity = next;
            }
            grid.branches[grid.nbranch++] = b;
        } else { status = GS_ERR_PARSE; break; }
    }
    if (ferror(file)) status = GS_ERR_IO;
    if (fclose(file) != 0 && status == GS_OK) status = GS_ERR_IO;
    if (status == GS_OK && !have_base) { grid.base_mva = 100.0; }
    if (status == GS_OK) status = gs_case_validate(&grid);
    if (status != GS_OK) {
        (void)line_number;
        gs_case_free(&grid);
        return status;
    }
    *out = grid;
    return GS_OK;
}