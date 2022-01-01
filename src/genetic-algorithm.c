#include "genetic-algorithm.h"

static Individual select_individual_with_replacement(const Individual *individuals, const unsigned n_individuals);
static Individual individual_crossover(const Individual g1, const Individual g2);
static Individual mutate_individual(const Individual g);
