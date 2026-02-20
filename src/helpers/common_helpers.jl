const GAM_INDEX = 1
const NS_INDEX = 2
const NR_INDEX = 3
const NE_INDEX = 4
const PASS_INDEX = 5

total_population(u) = u[NS_INDEX] + u[NR_INDEX] + u[NE_INDEX]
logistic_factor(N::Real, Cc::Real) = 1 - (N / Cc)
logistic_factor(u::AbstractVector, Cc::Real) = logistic_factor(total_population(u), Cc)

"""Return the population slice (sensitive, resistant, escape)."""
pop_fun(x) = x[NS_INDEX:NE_INDEX]

"""Return the passage counter from a state vector."""
pass_fun(x) = x[PASS_INDEX]

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
