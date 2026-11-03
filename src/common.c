#include "gridsolve.h"

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

void gs_log(int level, int configured_level, const char *format, ...)
{
    static const char *names[] = {"ERROR", "WARN", "INFO", "DEBUG"};
    va_list args;
    if (level > configured_level || level < GS_ERROR || level > GS_DEBUG) return;
    (void)fprintf(stderr, "[%s] ", names[level]);
    va_start(args, format);
    (void)vfprintf(stderr, format, args);
    va_end(args);
    (void)fputc('\n', stderr);
}

const char *gs_status_string(GsStatus status)
{
    switch (status) {
    case GS_OK: return "success";
    case GS_ERR_IO: return "I/O error";
    case GS_ERR_PARSE: return "parse error";
    case GS_ERR_INVALID: return "invalid case or argument";
    case GS_ERR_NOMEM: return "out of memory";
    case GS_ERR_SINGULAR: return "singular or ill-conditioned system";
    case GS_ERR_NOCONV: return "power flow did not converge";
    case GS_ERR_ISLAND: return "network is islanded";
    }
    return "unknown error";
}

GsOptions gs_default_options(void)
{
    GsOptions options = {1e-8, 20, 0, 0, GS_WARN};
    return options;
}

void gs_result_free(GsResult *result)
{
    if (result == NULL) return;
    free(result->vm);
    free(result->va_rad);
    free(result->p_inj);
    free(result->q_inj);
    free(result->mismatch_history);
    *result = (GsResult){0};
}