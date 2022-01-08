#include "RKF78.h"
#include "equations.h"
#include <float.h>
#include <math.h>

/* Elliot sigmoid Θ-scaled, σ-strengthened, and δ-displaced.
 */
static double sigmoid(const double x, const double sigma, const double delta);

/* Directed Elliot sigmoid.
 */
static double sigmoid_dir(const double x, const double mu, const double sigma, const double delta);

/* Dispersal by social coping.
 *
 * The function is designed so that the dispersal response of the population of
 * birds generically increases when the population numbers at the patch
 * diminish.
 */
static double model_dispersal(const double x, const double mu, const double sigma, const double delta);

/* Model to fit to the second epoch to check the hypothesis that the migration
 * of Andouins occurs with social copying.
 */
double model_equation(const double x, const Phenotype *const p);

/* Adaptation of model_equation function to fit the signature required by RKF78.
 */
void model_ode(double t, double x, double *result, void *p);

/* Intrinsic growth rate over the carrying capacity (birds^2/year).
 *
 * Estimated with the first epoch data to be 0.000024382635446.
 *
 * Initially defined as β = γ/K, where K = 16651.2696 is the carrying capacity.
 */
static const double beta = 0.000024382635446;

/* Scale of the sigmoid function.
 *
 * It is related with the order of magnitude of the carrying capacity K.
 */
static const double theta = 1000.0;

static double sigmoid(const double x, const double sigma, const double delta) {
    const double numerator = sigma * (x - delta);
    const double denominator = theta + sigma * fabs(x - delta);

    return numerator / denominator;
}

static double sigmoid_dir(const double x, const double mu, const double sigma, const double delta) {
    const double dir = mu * ((theta + sigma * delta) / (2 * theta + sigma * delta)) * (1 - x / delta) + x / delta;

    return dir * sigmoid(x, sigma, delta);
}

static double model_dispersal(const double x, const double mu, const double sigma, const double delta) {
    const double numerator = 1 - (x > delta ? sigmoid(x, sigma, delta) : sigmoid_dir(x, mu, sigma, delta));
    const double denominator = 1 - sigmoid_dir(0, mu, sigma, delta);

    return numerator / denominator;
}

double model_equation(const double x, const Phenotype *const p) {
    const double base = (p->phi * x) - (beta * x * x);

    return base - p->lambda * model_dispersal(x, p->mu, p->sigma, p->delta);
}

void model_ode(double __attribute__((unused)) t, double x, double *result, void *p) {
    *result = model_equation(x, p);
}

int model_prediction(const double x0, double *const x, const unsigned length, const Phenotype *const p) {
    double t = 0.0;
    double y = x0;
    double step = 1.0e-2;
    double error;
    const double step_min = 1.0e-3;
    const double step_max = 1.0e-2;
    const double tolerance = 1.0e-8;

    // Variables iter and t_end both count the same thing, but currently iter
    // does so as an unsigned integer to index the x vector, while t_end is a
    // double to easily compare it to t + step.
    //
    // This is done to avoid having to convert from integer to double or back.
    //
    // TODO: Check if it is worth it, performance wise.
    unsigned iter = 1;
    for(double t_end = 1; t_end < length; t_end++) {
        while(t + step < t_end) {
            int result = RKF78(&t, &y, &step, &error, step_min, step_max, tolerance, (void *) p, model_ode);
            if(result != 0) {
                return result;
            }
            if(!isnormal(y)) {
                return 1;
            }
        }
        step = t_end - t;

        int result = RKF78(&t, &y, &step, &error, step_min, step_max, tolerance, (void *) p, model_ode);
        if(result != 0) {
            return result;
        }
        if(!isnormal(y)) {
            return 1;
        }

        x[iter] = y;
        iter++;
    }
    x[0] = x0;

    return 0;
}

double get_phenotype_fitness(const Phenotype p) {
    double x[12] = { 15329.0 };
    int err = model_prediction(x[0], x, 12, &p);
    if(err != 0) {
        return DBL_MAX;
    }

    double fitness = 0.0;
    const double y[12] = { 15329.0, 14177.0, 13031.0, 9762.0, 11271.0, 8688.0, 7571.0, 6983.0, 4778.0, 2067.0, 1586.0, 793.0 };
    const double w[12] = {     1.0,     1.0,     1.0,    0.0,     1.0,    1.0,    1.0,    1.0,    3.0,    3.0,    3.0,  10.0 };
    for(unsigned char iter = 1; iter < 12; iter++) {
        const double tmp = w[iter] * (y[iter] - x[iter]) * (y[iter] - x[iter]);
        fitness += tmp;
        // if(tmp > fitness) {
            // fitness = tmp;
        // }
    }

    return fitness;
}
