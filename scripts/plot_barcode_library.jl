# Run with a Julia environment containing BarcodeModels, DataFrames, and RCall,
# and an R installation with ggplot2. From the repository root:
# julia --project=../julia_environments/bc_pkg_env_1.12.5 scripts/plot_barcode_library.jl 1000 output_directory

using BarcodeModels
using DataFrames
using Random
using RCall

# n0 is the number of starting cells. Several cells can receive the same barcode.
n0 = length(ARGS) >= 1 ? parse(Int64, ARGS[1]) : 1_000
output_dir = length(ARGS) >= 2 ? abspath(ARGS[2]) : pwd()
bc_unif_values = [0.0, 0.01, 0.1, 10.0]
Nbc_values = Int64[100, 1_000, 10_000]
seed = 42

n0 > 0 || throw(ArgumentError("n0 must be positive"))
Random.seed!(seed)

# Each parameter combination uses one fresh library and one sample of n0 cells.
barcode_counts = DataFrame(Nbc = Int64[], bc_unif = Float64[], cells_per_barcode = Int64[])
barcode_summary = DataFrame(Nbc = Int64[], bc_unif = Float64[], distinct_barcodes = Int64[])

for Nbc in Nbc_values, bc_unif in bc_unif_values
    sampled_barcodes = BarcodeModels._sample_initial_barcodes(n0;
        skew_lib = true, bc_unif = bc_unif, Nbc = Nbc)

    counts = zeros(Int64, Nbc)
    for barcode in sampled_barcodes
        counts[barcode] += 1
    end

    observed_counts = filter(>(0), counts)
    append!(barcode_counts, DataFrame(
        Nbc = fill(Nbc, length(observed_counts)),
        bc_unif = fill(bc_unif, length(observed_counts)),
        cells_per_barcode = observed_counts,
    ))
    push!(barcode_summary, (Nbc, bc_unif, length(observed_counts)))
end

mkpath(output_dir)
@rput barcode_counts barcode_summary n0 seed output_dir

R"""
library(ggplot2)

library_sizes <- sort(unique(barcode_summary$Nbc))
uniformity_values <- sort(unique(barcode_summary$bc_unif))
library_labels <- paste0("Nbc = ", format(library_sizes, big.mark = ",", trim = TRUE))
uniformity_labels <- paste0("bc_unif = ", uniformity_values)
barcode_counts$library_size <- factor(barcode_counts$Nbc, library_sizes, library_labels)
barcode_counts$uniformity <- factor(barcode_counts$bc_unif, uniformity_values, uniformity_labels)
barcode_summary$library_size <- factor(barcode_summary$Nbc, library_sizes, library_labels)
barcode_summary$uniformity <- factor(barcode_summary$bc_unif, uniformity_values)

plot_subtitle <- paste0(format(n0, big.mark = ","), " starting cells per setting  |  Seed ", seed)
plot_theme <- theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), plot.title.position = "plot")

barcode_histogram <- ggplot(barcode_counts, aes(x = cells_per_barcode)) +
    geom_histogram(binwidth = 1, boundary = 0, closed = "left", fill = "#247687") +
    scale_x_continuous(trans = "log2", breaks = 2^seq(0, ceiling(log2(n0)), by = 2)) +
    facet_wrap(vars(library_size, uniformity), ncol = length(uniformity_values), scales = "free_y") +
    labs(title = "Cells per observed barcode", subtitle = plot_subtitle,
         x = "Cells per barcode (log2 scale)", y = "Number of barcodes",
         caption = "Bins cover 1, 2 to 3, 4 to 7 cells and so on. The y-axis varies by panel.") +
    plot_theme

distinct_bar_plot <- ggplot(barcode_summary, aes(x = uniformity, y = distinct_barcodes)) +
    geom_col(fill = "#247687", width = 0.7) +
    geom_text(aes(label = distinct_barcodes), vjust = -0.4, size = 3.4) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.13))) +
    facet_wrap(vars(library_size), nrow = 1) +
    labs(title = "Distinct barcodes sampled", subtitle = plot_subtitle,
         x = "bc_unif", y = "Number of distinct barcodes") +
    plot_theme

ggsave(file.path(output_dir, "barcode_count_histograms.png"), barcode_histogram,
       width = 11, height = 7.5, dpi = 180)
ggsave(file.path(output_dir, "distinct_barcodes.png"), distinct_bar_plot,
       width = 9, height = 4, dpi = 180)
"""

println("Saved both figures to ", output_dir)
