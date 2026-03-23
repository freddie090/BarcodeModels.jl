const RESPOP_GAM_INDEX = 1
const RESPOP_NS_INDEX = 2
const RESPOP_NR_INDEX = 3
const RESPOP_NE_INDEX = 4
const RESPOP_PASS_INDEX = 5

const RESDMG_GAM_INDEX = 1
const RESDMG_NS_INDEX = 2
const RESDMG_NDS_INDEX = 3
const RESDMG_NDR_INDEX = 4
const RESDMG_NR_INDEX = 5
const RESDMG_PASS_INDEX = 6

respop_total_population(u) = u[RESPOP_NS_INDEX] + u[RESPOP_NR_INDEX] + u[RESPOP_NE_INDEX]
resdmg_total_population(u) = u[RESDMG_NS_INDEX] + u[RESDMG_NDS_INDEX] + u[RESDMG_NDR_INDEX] + u[RESDMG_NR_INDEX]

logistic_factor(N::Real, Cc::Real) = 1 - (N / Cc)
respop_logistic_factor(u::AbstractVector, Cc::Real) = logistic_factor(respop_total_population(u), Cc)
resdmg_logistic_factor(u::AbstractVector, Cc::Real) = logistic_factor(resdmg_total_population(u), Cc)

"""Return the ResPop population slice (sensitive, resistant, escape)."""
respop_pop_fun(x) = x[RESPOP_NS_INDEX:RESPOP_NE_INDEX]

"""Return the ResPop passage counter from a state vector."""
respop_pass_fun(x) = x[RESPOP_PASS_INDEX]

"""Return the ResDmg population slice (sensitive, sensitive-damaged, resistant-damaged, resistant)."""
resdmg_pop_fun(x) = x[RESDMG_NS_INDEX:RESDMG_NR_INDEX]

"""Return the ResDmg passage counter from a state vector."""
resdmg_pass_fun(x) = x[RESDMG_PASS_INDEX]

"""Draw multivariate hypergeometric counts (uses RNG)."""
function multivariate_hypergeometric_draw(population_sizes, num_samples)
    remaining_samples = num_samples
    total_pop = sum(population_sizes)
    draws = zeros(Int, length(population_sizes))

    for i in 1:(length(population_sizes) - 1)
        if remaining_samples > 0
            draws[i] = rand(Hypergeometric(population_sizes[i],
                                            total_pop - population_sizes[i],
                                            remaining_samples))
            remaining_samples -= draws[i]
            total_pop -= population_sizes[i]
        end
    end

    draws[end] = remaining_samples
    return draws
end

"""Return the index of `xvec` closest to `x` (ties pick the lower index)."""
function find_closest(x, xvec)
    n = length(xvec)
    idx = searchsortedfirst(xvec, x)

    if idx == 1
        return 1
    elseif idx == n
        return n
    else
        dist_prev = x - xvec[idx - 1]
        dist_next = xvec[idx] - x
        return dist_prev < dist_next ? idx - 1 : idx
    end
end
