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

const RESPOP_INVIVO_GAM_INDEX = 1
const RESPOP_INVIVO_NS_EG0_INDEX = 2
const RESPOP_INVIVO_NS_EG1_INDEX = 3
const RESPOP_INVIVO_NR_EG0_INDEX = 4
const RESPOP_INVIVO_NR_EG1_INDEX = 5
const RESPOP_INVIVO_NE_EG0_INDEX = 6
const RESPOP_INVIVO_NE_EG1_INDEX = 7
const RESPOP_INVIVO_PASS_INDEX = 8

respop_invivo_total_population(u) = u[2] + u[3] + u[4] + u[5] + u[6] + u[7]

function _validate_engraftment_probabilities(pEG, sEG)
    isfinite(pEG) && 0.0 <= pEG <= 1.0 || error("pEG must be finite and between 0 and 1.")
    isfinite(sEG) && 0.0 <= sEG <= 1.0 || error("sEG must be finite and between 0 and 1.")
    return nothing
end

"""Thin living EG-tagged cells independently; return survivors and realized counts."""
function engraftment_selection(cells::AbstractVector, pEG::Real, sEG::Real)
    _validate_engraftment_probabilities(pEG, sEG)
    keep = [cell.alive && rand() < (cell.EG ? pEG : pEG * (1-sEG)) for cell in cells]
    selected = cells[keep]
    nEG1 = count(cell -> cell.EG, selected)
    return selected, Dict("N_engraft" => length(selected),
        "nEG0_engraft" => length(selected)-nEG1, "nEG1_engraft" => nEG1)
end

"""Binomial engraftment thinning of nonnegative integer-like EG-stratified counts."""
function engraftment_selection(nEG0::AbstractVector, nEG1::AbstractVector, pEG::Real, sEG::Real)
    _validate_engraftment_probabilities(pEG, sEG)
    length(nEG0) == length(nEG1) || error("nEG0 and nEG1 must have matching lengths.")
    all(x -> isfinite(x) && x >= 0, nEG0) || error("nEG0 entries must be finite and nonnegative.")
    all(x -> isfinite(x) && x >= 0, nEG1) || error("nEG1 entries must be finite and nonnegative.")
    eg0 = [rand(Binomial(round(Int, n), pEG*(1-sEG))) for n in nEG0]
    eg1 = [rand(Binomial(round(Int, n), pEG)) for n in nEG1]
    return eg0, eg1, Dict("N_engraft" => sum(eg0)+sum(eg1),
        "nEG0_engraft" => sum(eg0), "nEG1_engraft" => sum(eg1))
end
