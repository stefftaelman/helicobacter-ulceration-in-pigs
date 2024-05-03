
using CSV, Colors, DataFrames, Plots

COLOR_SCHEME = [
    colorant"#00A98F", # biolizard green
    colorant"#FFC000", # gold
    colorant"#1565A9", # blue
    colorant"#C00000", # red
    colorant"#0D0D0D"  # black
]

# set paths
pre_qc = "./data/mqc_phred_scores_pre.tsv"
post_qc = "./data/mqc_phred_scores_post.tsv"

# load in tsv files
pre_qc_df = CSV.read(pre_qc, DataFrame, drop=[:Sample])
post_qc_df = CSV.read(post_qc, DataFrame, drop=[:Sample])

# make line plots plotting each sample's phred scores at the positions in their column headers

plot()
for i in 1:size(pre_qc_df, 1)
    phred_scores=Vector(pre_qc_df[i, :])
    c = any(phred_scores .< 20) ? COLOR_SCHEME[4] : any(phred_scores .< 23) ? COLOR_SCHEME[2] : COLOR_SCHEME[1]
    plot!(parse.(Int, names(pre_qc_df)), phred_scores, label="", color=c)
    current()
end
display(plot!())
ylims!(0, 39)
yticks!(0:5:40)
savefig("./figures/phred_scrores_preQC.svg")

plot()
for i in 1:size(post_qc_df, 1)
    phred_scores=Vector(post_qc_df[i, :])
    c = any(phred_scores .< 20) ? COLOR_SCHEME[4] : any(phred_scores .< 23) ? COLOR_SCHEME[2] : COLOR_SCHEME[1]
    plot!(parse.(Int, names(post_qc_df)), phred_scores, label="", color=c)
    current()
end
display(plot!())
ylims!(0, 39)
yticks!(0:5:40)
savefig("./figures/phred_scrores_postQC.svg")