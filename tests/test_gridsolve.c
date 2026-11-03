#include "gridsolve.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int tests_run;
static int tests_failed;
static double hand_max_voltage_error;
static double hand_max_angle_error;
static double max_branch_p_error;
static double max_loss_error;
static double ieee_max_voltage_error;
static double ieee_max_angle_error;

#define CHECK(expr) do { ++tests_run; if (!(expr)) { \
    (void)fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #expr); ++tests_failed; \
} } while (0)

static void two_terminal_hand_solution(double p, double q, double x,
                                       double *vm, double *angle)
{
    double a = 1.0 - 2.0*q*x;
    double discriminant = a*a - 4.0*x*x*(p*p+q*q);
    double y = (a + sqrt(discriminant))/2.0;
    *vm = sqrt(y);
    *angle = atan2(-p*x, y+q*x);
}

static void test_lu(void)
{
    double complex a[] = {2,1,1,3};
    double complex b[] = {1,2};
    double complex singular[] = {1,2,2,4};
    double complex c[] = {1,1};
    CHECK(gs_lu_solve(2, a, b) == GS_OK);
    CHECK(fabs(creal(b[0])-0.2) < 1e-14 && fabs(creal(b[1])-0.6) < 1e-14);
    CHECK(gs_lu_solve(2, singular, c) == GS_ERR_SINGULAR);
}

static void test_ybus_and_two_bus(void)
{
    GsCase grid = {0};
    GsOptions options = gs_default_options();
    GsResult result = {0};
    double complex *y = NULL;
    double complex sentinel;
    double expected_vm, expected_angle;
    size_t i, j;
    y = &sentinel;
    CHECK(gs_ybus(NULL, &y) == GS_ERR_INVALID);
    CHECK(y == NULL);
    CHECK(gs_case_load("data/two_bus.grid", &grid) == GS_OK);
    CHECK(gs_ybus(&grid, &y) == GS_OK);
    if (y != NULL) {
        for (i = 0; i < grid.nbus; ++i) {
            double complex sum = 0;
            for (j = 0; j < grid.nbus; ++j) {
                CHECK(cabs(y[i*grid.nbus+j]-y[j*grid.nbus+i]) < 1e-13);
                sum += y[i*grid.nbus+j];
            }
            CHECK(cabs(sum) < 1e-13);
        }
    }
    CHECK(gs_solve(&grid, &options, &result) == GS_OK);
    if (result.vm != NULL) {
        /* From S=V*conj((1-V)/(jX)), y=|V|^2 solves
           y^2+(2QX-1)y+X^2(P^2+Q^2)=0; use the high-voltage root. */
        two_terminal_hand_solution(0.5, 0.2, 0.1, &expected_vm, &expected_angle);
        hand_max_voltage_error = fmax(hand_max_voltage_error, fabs(result.vm[1]-expected_vm));
        hand_max_angle_error = fmax(hand_max_angle_error, fabs(result.va_rad[1]-expected_angle)*180.0/acos(-1.0));
        CHECK(fabs(result.vm[1]-expected_vm) < 1e-8);
        CHECK(fabs(result.va_rad[1]-expected_angle) < 1e-8);
        CHECK(fabs(result.total_loss_p) < 1e-7);
        CHECK(result.total_loss_q > 0.0);
        {
            GsFlow flow;
            double q_loss = 0.1*(0.5*0.5+0.2*0.2)/(result.vm[1]*result.vm[1])*100.0;
            CHECK(gs_branch_flows(&grid, &result, &flow) == GS_OK);
            CHECK(fabs(flow.p_from-50.0) < 1e-7);
            CHECK(fabs(flow.q_from-(20.0+q_loss)) < 1e-7);
            CHECK(fabs(flow.p_to+50.0) < 1e-7);
            CHECK(fabs(flow.q_to+20.0) < 1e-7);
        }
    }
    free(y); gs_result_free(&result); gs_case_free(&grid);
}

static void test_three_bus(void)
{
    GsCase grid = {0};
    GsResult result = {0};
    double expected_vm, expected_angle;
    CHECK(gs_case_load("data/three_bus.grid", &grid) == GS_OK);
    CHECK(gs_solve(&grid, NULL, &result) == GS_OK);
    if (result.vm != NULL) {
        two_terminal_hand_solution(0.3, 0.1, 0.2, &expected_vm, &expected_angle);
        hand_max_voltage_error = fmax(hand_max_voltage_error, fabs(result.vm[2]-expected_vm));
        hand_max_angle_error = fmax(hand_max_angle_error, fabs(result.va_rad[2]-expected_angle)*180.0/acos(-1.0));
        CHECK(fabs(result.vm[2]-expected_vm) < 1e-8);
        CHECK(fabs(result.va_rad[2]-expected_angle) < 1e-8);
        /* Equal series currents make the middle-bus phasor the endpoint midpoint. */
        {
            double complex v1 = 1.0;
            double complex v3 = result.vm[2]*(cos(result.va_rad[2])+I*sin(result.va_rad[2]));
            double complex v2 = result.vm[1]*(cos(result.va_rad[1])+I*sin(result.va_rad[1]));
            CHECK(cabs(v2-(v1+v3)/2.0) < 1e-8);
        }
    }
    gs_result_free(&result); gs_case_free(&grid);
}

static void test_q_limits_and_quadratic_rate(void)
{
    GsCase grid = {0};
    GsResult result = {0};
    GsOptions options = gs_default_options();
    size_t i;
    int quadratic = 0;
    CHECK(gs_case_load("data/three_bus.grid", &grid) == GS_OK);
    if (grid.buses != NULL) {
        grid.buses[1].type = GS_PV;
        grid.buses[1].qmin = 0.0;
        grid.buses[1].qmax = 0.0;
        options.enforce_q_limits = 1;
        CHECK(gs_solve(&grid, &options, &result) == GS_OK);
        if (result.vm != NULL) {
            CHECK(fabs(result.q_inj[1]*grid.base_mva) < 1e-7);
            gs_result_free(&result);
        }
        options.enforce_q_limits = 0;
        options.tolerance = 1e-13;
        CHECK(gs_solve(&grid, &options, &result) == GS_OK);
        if (result.vm != NULL) CHECK(fabs(result.vm[1]-grid.buses[1].vm) < 1e-10);
        if (result.mismatch_history != NULL) {
            for (i = 1; i < result.mismatch_count; ++i) {
                double previous = result.mismatch_history[i-1];
                double current = result.mismatch_history[i];
                if (previous > 0.0 && previous < 1e-2 && current <= 100.0*previous*previous)
                    quadratic = 1;
            }
            CHECK(quadratic);
        } else CHECK(0);
    }
    gs_result_free(&result); gs_case_free(&grid);
}

static void test_ieee14_web_reference(void)
{
    /* Reference generated by the browser solver at the URL documented in README.md. */
    static const double reference_vm[] = {
        1.06, 1.045, 1.01, 1.0176708575976607, 1.0195138638602699,
        1.07, 1.061519537688715, 1.09, 1.0559317263657604,
        1.0509846305368538, 1.0569065225314898, 1.0551885642273768,
        1.0503817153543678, 1.0355299506967024
    };
    static const double reference_angle_deg[] = {
        0.0, -4.98258900932882, -12.72509967768685, -10.3129008692666,
        -8.77385367338953, -14.220945956274239, -13.35962717444426,
        -13.35962717444426, -14.938521236574106, -15.097288304257562,
        -14.790621618351427, -15.075584017189437, -15.15627586261134,
        -16.03364430874026
    };
    static const double reference_branch_p_from[] = {
        156.8828865215467, 75.5103798261868, 73.23757808578245,
        56.13149456708246, 41.51621353792567, -23.285690627479383,
        -61.158231082080846, 28.07417645303373, 16.079758302322414,
        44.087318760165026, 7.3532769320991065, 7.786066878767373,
        17.74797671657451, 0.0, 28.07417855701637, 5.22755058075175,
        9.426380300430015, -3.785324270061219, 1.6142577199079522,
        5.643851755753526
    };
    static const double reference_total_loss_mw = 13.393271730065992;
    GsCase grid = {0};
    GsResult result = {0};
    GsFlow flows[20];
    size_t i;
    CHECK(gs_case_load("data/case14.grid", &grid) == GS_OK);
    CHECK(grid.nbus == 14 && grid.nbranch == 20);
    CHECK(gs_solve(&grid, NULL, &result) == GS_OK);
    if (result.vm != NULL) {
        CHECK(gs_branch_flows(&grid, &result, flows) == GS_OK);
        for (i = 0; i < 14; ++i) {
            double vm_error = fabs(result.vm[i]-reference_vm[i]);
            double angle_error = fabs(result.va_rad[i]*180.0/acos(-1.0)-reference_angle_deg[i]);
            if (vm_error > ieee_max_voltage_error) ieee_max_voltage_error = vm_error;
            if (angle_error > ieee_max_angle_error) ieee_max_angle_error = angle_error;
            CHECK(vm_error < 1e-6);
            CHECK(angle_error < 1e-4);
        }
        for (i = 0; i < 20; ++i) {
            double error = fabs(flows[i].p_from-reference_branch_p_from[i]);
            if (error > max_branch_p_error) max_branch_p_error = error;
            CHECK(error < 1e-4);
        }
        max_loss_error = fabs(result.total_loss_p-reference_total_loss_mw);
        CHECK(max_loss_error < 1e-4);
    }
    gs_result_free(&result);
    gs_case_free(&grid);
}

static void test_failures_and_balance(void)
{
    GsCase grid = {0}, broken = {0};
    GsResult result = {0};
    double p_balance = 0.0, q_balance = 0.0;
    size_t i;
    CHECK(gs_case_load("data/two_bus.grid", &grid) == GS_OK);
    CHECK(gs_solve(&grid, NULL, &result) == GS_OK);
    if (result.vm != NULL) {
        for (i = 0; i < grid.nbus; ++i) {
            p_balance += result.p_inj[i]*grid.base_mva;
            q_balance += result.q_inj[i]*grid.base_mva;
        }
        CHECK(fabs(p_balance-result.total_loss_p) < 1e-6);
        CHECK(fabs(q_balance-result.total_loss_q) < 1e-6);
    }
    gs_result_free(&result);
    CHECK(gs_case_clone(&grid, &broken) == GS_OK);
    if (broken.branches != NULL) {
        broken.branches[0].status = 0;
        CHECK(gs_case_validate(&broken) == GS_ERR_ISLAND);
        broken.branches[0].status = 1;
        broken.branches[0].x = -0.1;
        CHECK(gs_case_validate(&broken) == GS_ERR_INVALID);
    }
    gs_case_free(&broken);
    CHECK(gs_case_clone(&grid, &broken) == GS_OK);
    if (broken.buses != NULL) {
        broken.buses[1].pd = 100000.0;
        broken.buses[1].qd = 40000.0;
        CHECK(gs_solve(&broken, NULL, &result) == GS_ERR_NOCONV);
    }
    gs_result_free(&result); gs_case_free(&broken); gs_case_free(&grid);
    CHECK(gs_case_load("data/two_bus.grid", &grid) == GS_OK);
    if (grid.buses != NULL) {
        grid.buses[1].qd = -3.9;
        CHECK(gs_case_validate(&grid) == GS_OK);
    }
    gs_case_free(&grid);
}

static void test_contingency_and_sweep(void)
{
    GsCase grid = {0}, triangle = {0}, limited = {0};
    GsContingencyResult *events = NULL;
    size_t count = 0, i;
    double first, collapse;
    GsStatus status;
    CHECK(gs_case_load("data/three_bus.grid", &grid) == GS_OK);
    CHECK(gs_case_clone(&grid, &triangle) == GS_OK);
    if (triangle.branches != NULL) {
        GsBranch *grown = realloc(triangle.branches, 3*sizeof(*grown));
        CHECK(grown != NULL);
        if (grown != NULL) {
            triangle.branches = grown;
            triangle.branches[2] = (GsBranch){1,3,1,0.0,0.2,0.0,1.0,0.0,10.0};
            triangle.nbranch = 3;
            CHECK(gs_n1(&triangle, NULL, &events, &count) == GS_OK);
            CHECK(count == 3);
            if (events != NULL) {
                int saw_overload = 0, saw_connected = 0;
                for (i = 0; i < count; ++i) {
                    if (events[i].status == GS_OK) {
                        ++saw_connected;
                        if (events[i].overloads > 0) saw_overload = 1;
                    }
                    CHECK(events[i].branch_index >= 0 && events[i].branch_index < 3);
                }
                CHECK(saw_connected == 3);
                CHECK(saw_overload);
            }
            free(events); events = NULL; count = 0;
        }
    }
    CHECK(gs_case_load("data/two_bus.grid", &limited) == GS_OK);
    if (limited.buses != NULL) {
        limited.buses[1].vmin = 0.999;
        status = gs_sweep(&limited, NULL, 0.5, 0.6, 0.1, &first, &collapse);
        CHECK(status == GS_OK);
        CHECK(fabs(first-0.5) < 1e-12);
        CHECK(isnan(collapse));
        limited.buses[1].vmin = 0.1;
        limited.buses[1].pd = 1000.0;
        limited.buses[1].qd = 400.0;
        status = gs_sweep(&limited, NULL, 0.5, 1.5, 0.5, &first, &collapse);
        CHECK(status == GS_OK);
        CHECK(isfinite(collapse) && collapse >= 0.5 && collapse <= 1.5);
    }
    gs_case_free(&grid); gs_case_free(&triangle); gs_case_free(&limited);
}

int main(void)
{
    test_lu();
    test_ybus_and_two_bus();
    test_three_bus();
    test_q_limits_and_quadratic_rate();
    test_ieee14_web_reference();
    test_failures_and_balance();
    test_contingency_and_sweep();
    (void)printf("%d checks, %d failures\n", tests_run, tests_failed);
    (void)printf("hand-reference max errors: |V| %.12e pu, angle %.12e deg\n",
                 hand_max_voltage_error, hand_max_angle_error);
    (void)printf("IEEE14 browser-reference max errors: |V| %.12e pu, angle %.12e deg, branch P %.12e MW, total loss %.12e MW\n",
                 ieee_max_voltage_error, ieee_max_angle_error, max_branch_p_error, max_loss_error);
    return tests_failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}