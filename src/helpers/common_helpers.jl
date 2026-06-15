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

const RESPOP_INVIVO_GAM_INDEX = 1
const RESPOP_INVIVO_NS_EG0_INDEX = 2
const RESPOP_INVIVO_NS_EG1_INDEX = 3
const RESPOP_INVIVO_NR_EG0_INDEX = 4
const RESPOP_INVIVO_NR_EG1_INDEX = 5
const RESPOP_INVIVO_NE_EG0_INDEX = 6
const RESPOP_INVIVO_NE_EG1_INDEX = 7
const RESPOP_INVIVO_PASS_INDEX = 8

respop_total_population(u) = u[RESPOP_NS_INDEX] + u[RESPOP_NR_INDEX] + u[RESPOP_NE_INDEX]
resdmg_total_population(u) = u[RESDMG_NS_INDEX] + u[RESDMG_NDS_INDEX] + u[RESDMG_NDR_INDEX] + u[RESDMG_NR_INDEX]
respop_invivo_total_population(u) =
    u[RESPOP_INVIVO_NS_EG0_INDEX] + u[RESPOP_INVIVO_NS_EG1_INDEX] +
    u[RESPOP_INVIVO_NR_EG0_INDEX] + u[RESPOP_INVIVO_NR_EG1_INDEX] +
    u[RESPOP_INVIVO_NE_EG0_INDEX] + u[RESPOP_INVIVO_NE_EG1_INDEX]

solver_suppression_enabled() =
    lowercase(get(ENV, "BARCODEMODELS_SUPPRESS_SOLVER", "true")) in ("1", "true", "yes", "on")

macro maybe_suppress_solver(expr)
    return esc(quote
        if solver_suppression_enabled()
            @suppress $expr
        else
            $expr
        end
    end)
end

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

"""
Apply in vivo engraftment selection to EG-tagged cells.

Returns `(selected_cells, stats)` where `stats` contains:
`N_engraft`, `nEG0_engraft`, `nEG1_engraft`.
"""
function engraftment_selection(cells::AbstractVector, pEG::Float64, sEG::Float64)
    n = length(cells)
    keep_mask = falses(n)
    nEG0_engraft = 0
    nEG1_engraft = 0
    p_eg0 = pEG * (1.0 - sEG)

    for i in eachindex(cells)
        p_keep = cells[i].EG ? pEG : p_eg0
        if rand() <= p_keep
            keep_mask[i] = true
            if cells[i].EG
                nEG1_engraft += 1
            else
                nEG0_engraft += 1
            end
        end
    end

    selected = cells[keep_mask]
    stats = Dict(
        "N_engraft" => nEG0_engraft + nEG1_engraft,
        "nEG0_engraft" => nEG0_engraft,
        "nEG1_engraft" => nEG1_engraft
    )
    return selected, stats
end

"""
Apply in vivo engraftment selection to EG0/EG1 count vectors.

`nEG0` and `nEG1` must be equal-length vectors of non-negative integer-like values.
Returns `(new_nEG0, new_nEG1, stats)` where counts are Int vectors and `stats`
contains `N_engraft`, `nEG0_engraft`, `nEG1_engraft`.
"""
function engraftment_selection(nEG0::AbstractVector, nEG1::AbstractVector, pEG::Float64, sEG::Float64)
    length(nEG0) == length(nEG1) || error("nEG0 and nEG1 must have matching lengths.")
    p_eg0 = pEG * (1.0 - sEG)

    new_nEG0 = Vector{Int64}(undef, length(nEG0))
    new_nEG1 = Vector{Int64}(undef, length(nEG1))

    for i in eachindex(nEG0)
        n0 = Int64(round(nEG0[i]))
        n1 = Int64(round(nEG1[i]))
        n0 >= 0 || error("nEG0 entries must be non-negative.")
        n1 >= 0 || error("nEG1 entries must be non-negative.")

        new_nEG0[i] = n0 == 0 ? 0 : rand(Binomial(n0, p_eg0))
        new_nEG1[i] = n1 == 0 ? 0 : rand(Binomial(n1, pEG))
    end

    nEG0_engraft = sum(new_nEG0)
    nEG1_engraft = sum(new_nEG1)
    stats = Dict(
        "N_engraft" => nEG0_engraft + nEG1_engraft,
        "nEG0_engraft" => nEG0_engraft,
        "nEG1_engraft" => nEG1_engraft
    )

    return new_nEG0, new_nEG1, stats
end
