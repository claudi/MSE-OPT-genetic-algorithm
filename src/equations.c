#include "equations.h"
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
