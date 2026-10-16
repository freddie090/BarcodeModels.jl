function model_meas_noise(u_outs::Array{Float64}, phi::Float64)

    if sum(u_outs .< 0) > 0
        noise_u_outs = u_outs
    else
        noise_u_outs = rand.(Normal.(u_outs, u_outs * phi))

        noise_u_outs[noise_u_outs .< 0] .= 0
        noise_u_outs = Int64.(round.(noise_u_outs, digits = -3))
    end

    return noise_u_outs
end
