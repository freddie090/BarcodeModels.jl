"""Select a behavior based on `drug_effect` (:d, :b, or :c)."""
function select_drug_effect(de, f_d, f_b, f_c)
    return de == :d ? f_d : de == :b ? f_b : de == :c ? f_c :
           error("Invalid drug_effect value: must be :b, :d, or :c")
end

safe_ratio(num, den) = abs(den) > eps(Float64) ? (num / den) : 0.0

function apply_drug_effect_death(b_ref, d_ref, b, d, Dc, gam, N, Cc, psi)
    lf = logistic_factor(N, Cc)
    b_mod = b * lf
    d_mod = (d + (safe_ratio(d, d_ref) * Dc * gam * (1 - psi))) * lf
    return b_mod, d_mod
end

function apply_drug_effect_birth(b_ref, d_ref, b, d, Dc, gam, N, Cc, psi)
    lf = logistic_factor(N, Cc)
    b_mod = (b - (safe_ratio(b, b_ref) * Dc * gam * (1 - psi))) * lf
    d_mod = d * lf
    return b_mod, d_mod
end

function apply_drug_effect_combined(b_ref, d_ref, b, d, Dc, gam, N, Cc, psi)
    lf = logistic_factor(N, Cc)
    b_eff = b - (safe_ratio(b, b_ref) * Dc * gam * (1 - psi))
    if b_eff >= 0
        b_mod = b_eff * lf
        d_mod = d * lf
    else
        b_mod = 0.0
        d_mod = (d + abs(b_eff)) * lf
    end
    return b_mod, d_mod
end
