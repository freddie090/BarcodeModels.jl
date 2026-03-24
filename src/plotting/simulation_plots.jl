using RCall

"""
    plot_simulation_outputs(
        hybrid_sol_df,
        hybrid_pop_df,
        treat_ons,
        treat_offs,
        plot_tmax,
        Nswitch;
        agent_based_sol_df = nothing,
        agent_based_pop_df = nothing,
        include_agent_based = true,
    )

Plot Hybrid-model simulation outputs, with an optional side-by-side Agent-Based panel.

When `include_agent_based=true`, both `agent_based_sol_df` and `agent_based_pop_df`
must be provided.
"""
function plot_simulation_outputs(hybrid_sol_df::DataFrame,
                                 hybrid_pop_df::DataFrame,
                                 treat_ons::Vector{Float64},
                                 treat_offs::Vector{Float64},
                                 plot_tmax::Float64,
                                 Nswitch::Int64;
                                 agent_based_sol_df::Union{Nothing, DataFrame}=nothing,
                                 agent_based_pop_df::Union{Nothing, DataFrame}=nothing,
                                 include_agent_based::Bool=true)

    if include_agent_based && (isnothing(agent_based_sol_df) || isnothing(agent_based_pop_df))
        error("When include_agent_based=true, both agent_based_sol_df and agent_based_pop_df must be provided.")
    end

    @rput hybrid_sol_df
    @rput hybrid_pop_df
    @rput treat_ons
    @rput treat_offs
    @rput plot_tmax
    @rput Nswitch
    @rput include_agent_based

    if include_agent_based
        @rput agent_based_sol_df
        @rput agent_based_pop_df
    else
        agent_based_sol_df = DataFrame()
        agent_based_pop_df = DataFrame()
        @rput agent_based_sol_df
        @rput agent_based_pop_df
    end

    R"""
        library(tidyr)
        library(dplyr)
        library(ggplot2)
        library(cowplot)

        theme_barcodemodels <- function(){
            theme_bw() +
                theme(
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_line(color = "grey90"),
                    legend.position = "bottom",
                    plot.title = element_text(face = "bold")
                )
        }

        make_treat_df_local <- function(treat_ons, treat_offs, tmax){
            if(length(treat_ons) == 0 || length(treat_offs) == 0){
                return(data.frame(start = numeric(0), end = numeric(0), treat = character(0)))
            }

            n <- min(length(treat_ons), length(treat_offs))
            starts <- pmax(treat_ons[seq_len(n)], 0.0)
            ends <- pmin(treat_offs[seq_len(n)], tmax)
            keep <- ends > starts

            data.frame(
                start = starts[keep],
                end = ends[keep],
                treat = "On"
            )
        }

        if("nDR" %in% colnames(hybrid_sol_df)){
            pop_pal <- c("black", "#B15940", "#6A3D9A", "#E31A1C", "#1F78B4")
        } else if("nRAB" %in% colnames(hybrid_sol_df)){
            pop_pal <- c("black", "#6A3D9A", "#E31A1C", "#33A02C", "#FF7F00", "#1F78B4")
        } else {
            pop_pal <- c("black", "#6A3D9A", "#E31A1C", "#1F78B4")
        }

        treat_df <- make_treat_df_local(treat_ons, treat_offs, plot_tmax)

        hybrid_sol_long <- hybrid_sol_df %>%
            tidyr::pivot_longer(cols = matches("^n[A-Z]+$"),
                                names_to = "pheno",
                                values_to = "n") %>%
            mutate(n = if_else(n < 1, 0, n))

        p_hybrid <- ggplot() +
            geom_rect(data = treat_df,
                      aes(ymin = 1.0, ymax = 1e+06,
                          xmin = start, xmax = end,
                          fill = treat), alpha = 0.1, colour = "black") +
            scale_fill_manual(values = c("On" = "red"), name = "Treatment") +
            geom_line(data = hybrid_sol_long,
                      aes(x = t, y = n, colour = pheno,
                          group = interaction(rep, pheno)), size = 1) +
            geom_line(data = hybrid_sol_long,
                      aes(x = t, y = N, colour = "N",
                          group = interaction(rep, pheno)),
                      linetype = "dashed", alpha = 0.2, size = 1.2) +
            geom_point(data = hybrid_pop_df, aes(x = t, y = N),
                       size = 3.0, alpha = 0.5, colour = "black") +
            geom_hline(yintercept = Nswitch, alpha = 0.5, linetype = "dotted") +
            scale_colour_manual(values = pop_pal) +
            scale_x_continuous(limits = c(0.0, plot_tmax)) +
            scale_y_log10() +
            ggtitle("Hybrid Model") +
            theme_barcodemodels()

        if("drug_pheno" %in% colnames(hybrid_sol_long)){
            p_hybrid <- p_hybrid + facet_wrap(~drug_pheno, ncol = 1)
        }

        if(!include_agent_based){
            print(p_hybrid)
        } else {
            agent_based_sol_long <- agent_based_sol_df %>%
                tidyr::pivot_longer(cols = matches("^n[A-Z]+$"),
                                    names_to = "pheno",
                                    values_to = "n") %>%
                mutate(n = if_else(n < 1, 0, n))

            p_agent_based <- ggplot() +
                geom_rect(data = treat_df,
                          aes(ymin = 1.0, ymax = 1e+06,
                              xmin = start, xmax = end,
                              fill = treat), alpha = 0.1, colour = "black") +
                scale_fill_manual(values = c("On" = "red"), name = "Treatment") +
                geom_line(data = agent_based_sol_long,
                          aes(x = t, y = n, colour = pheno,
                              group = interaction(rep, pheno)), size = 1) +
                geom_line(data = agent_based_sol_long,
                          aes(x = t, y = N, colour = "N",
                              group = interaction(rep, pheno)),
                          linetype = "dashed", alpha = 0.2, size = 1.2) +
                geom_point(data = agent_based_pop_df, aes(x = t, y = N),
                           size = 3.0, alpha = 0.5, colour = "black") +
                geom_hline(yintercept = Nswitch, alpha = 0.5, linetype = "dotted") +
                scale_colour_manual(values = pop_pal) +
                scale_x_continuous(limits = c(0.0, plot_tmax)) +
                scale_y_log10() +
                ggtitle("Agent-Based Model") +
                theme_barcodemodels()

            if("drug_pheno" %in% colnames(agent_based_sol_long)){
                p_agent_based <- p_agent_based + facet_wrap(~drug_pheno, ncol = 1)
            }

            print(cowplot::plot_grid(p_hybrid, p_agent_based, nrow = 1))
        }
    """

    return nothing
end
