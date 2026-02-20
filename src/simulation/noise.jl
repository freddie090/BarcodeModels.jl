function model_meas_noise(u_outs::Array{Float64}, phi::Float64,
                          run_colony::Bool = false)

    if sum(u_outs .< 0) > 0
        noise_u_outs = u_outs
    else
        if run_colony
            noise_u_outs = rand.(Normal.(u_outs[2:length(u_outs)],
                                         u_outs[2:length(u_outs)] * phi))
        else
            noise_u_outs = rand.(Normal.(u_outs, u_outs * phi))
        end

        noise_u_outs[noise_u_outs .< 0] .= 0
        noise_u_outs = Int64.(round.(noise_u_outs, digits = -3))

        if run_colony
            insert!(noise_u_outs, 1, u_outs[1])
        end
    end

    return noise_u_outs
end
