/* Contains the parameters for the model with equation
 *
 *     dx/dt = γx - βx^2 - λΨ(x, μ, σ, δ)
 *
 * where x(t) is the number of birds at time t.
 *
 * The function Ψ corresponds to the dispersal function.
 */
typedef struct {
    /* Intrinsic growth rate (birds^2/year).
     *
     * Using the first epoch it can be estimated to be γ = 0.406001835194.
     */
    double phi;
    /* Non-linear dispersal rate (birds/year).
     *
     * Must be lower than the neat population growth rate α, which can be estimated with the first epoch to be α = 0.3489494085776018.
     */
    double lambda;
    /* Determines the sign of the derivative of the dispersal function Ψ at x = 0.
     *
     *  - If μ is greater than 1, the sign will be negative,
     *  - If μ is lower than 1, the sign will be positive,
     *  - If μ is 1, the value of the derivative will be 0.
     *
     * Must be a positive number.
     */
    double mu;
    /* Determines the slopes of the sigmoids.
     * If σ ≈ 600, it approximates a Heaviside function.
     *
     * Must be a positive number.
     */
    double sigma;
    /* Point of change of concavity of Ψ(x, μ, σ, δ).
     *
     * Must be a positive number.
     */
    double delta;
} Phenotype;

/* Model to fit to the second epoch to check  the hypothesis that the migration
 * of Andouins occurs with social copying.
 */
double model_equation(const double x, const Phenotype *const p);
