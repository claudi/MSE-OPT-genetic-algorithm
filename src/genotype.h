#pragma once
#include "equations.h"

/* Structure containing the discretisation of the ODE parameters as unsigned
 * integers of appropriate length.
 *
 * This is done for compatibility with Holland's Convergence Theorem, which
 * works under the setting of genes consisting in unsigned integers expressed
 * in binary.
 *
 * For each parameter we need the theoretical range, an effective search range
 * and a reasonable precision, thus fixing the range and discretisation
 * formula for the genotype.
 *
 * The structure is packed as to occupy the least space possible, and in total
 * we achieve a `sizeof(Phenotype) = 16`, while `sizeof(Genotype) = 40`.
 */
typedef struct {
    /* - Theoretical search range: (-∞, α], where α = 0.3489494085776018 is the
     * neat population growth rate estimated from the first epoch.
     *
     * - Effective search range: [-100, 0.35].
     *
     * - Discretisation step: 1e-8.
     *
     * - Integer search range: [0, 2^34 - 1].
     */
    unsigned long phi : 34;
    /* - Theoretical search range: (0, ∞).
     *
     * - Effective search range: [0, 3000].
     *
     * - Discretisation step: 1e-4.
     *
     * - Integer search range: [0, 2^25 - 1].
     */
    unsigned lambda : 25;
    /* - Theoretical search range: (0, ∞).
     *
     * - Effective search range: [0, 20].
     *
     * - Discretisation step: 1e-6.
     *
     * - Integer search range: [0, 2^25 - 1].
     */
    unsigned mu : 25;
    /* - Theoretical search range: (0, ∞).
     *
     * - Effective search range: [0, 1000].
     *
     * - Discretisation step: 1e-2.
     *
     * - Integer search range: [0, 2^17 - 1].
     */
    unsigned sigma : 17;
    /* - Theoretical search range: (0, ∞).
     *
     * - Effective search range: [0, 25000].
     *
     * - Discretisation step: 1.
     *
     * - Integer search range: [0, 2^15 - 1].
     */
    unsigned short delta : 15;
} Genotype;

/* Function to convert from discrertised coefficients, used by the genetic
 * algorithm as unsigned integers, to floating point numbers.
 */
Phenotype genoype_to_phenotype(const Genotype g);

Genotype get_random_genotype();
Genotype genotype_crossover(const Genotype g1, const Genotype g2);
Genotype mutate_genotype(const Genotype g);
