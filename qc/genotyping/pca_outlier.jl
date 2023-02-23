#!/usr/bin/env julia
using CSV,CairoMakie,DataFrames,Statistics

pcafile = ARGS[1] # e.g. "/u/home/d/danieldu/shared-gandalm/GenomicDatasets/PsychENCODE/Genotypes/IndividualStudies/UCLA_ASD_synapse/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_EUR_brainCTP_ldprune_pca20.eigenvec"
outfile = ARGS[2] # e.g. "/u/home/d/danieldu/shared-gandalm/GenomicDatasets/PsychENCODE/Genotypes/IndividualStudies/UCLA_ASD_synapse/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_EUR_brainCTP_ldprune_pca20"

function getoutliers(pcadf::AbstractDataFrame, outlierstd::Real = 3)
    outliers = String[]
    for i in 1:3
        pc = pcadf[:, i + 2]
        pcoutlier = pcadf[abs.(pc .- mean(pc)) .> outlierstd * std(pc), 1]
        outliers = vcat(outliers, pcoutlier)
    end
    outliers = Set(outliers)
    outliers = [outliers...]
    df = DataFrame(:FamID => outliers,
                   :SampleID => outliers)
end

function graphpcs(pcadf::AbstractDataFrame, colind1::Integer, colind2::Integer, outlierstd::Real = 2)
    f = Figure()
    ax = Axis(f[1,1])
    colors = [:blue, :yellow, :orange, :red]
    pc1 = pcadf[:, colind1]
    pc2 = pcadf[:, colind2]
    meanpc1 = mean(pc1)
    meanpc2 = mean(pc2)
    stdpc1 = std(pc1)
    stdpc2 = std(pc2)
    for i in 0:3
        samplestoplot = pcadf[abs.(pc1 .- meanpc1) .> i * stdpc1 .||
                              abs.(pc2 .- meanpc2) .> i * stdpc2, [colind1, colind2]]
        scatter!(ax, samplestoplot[:, 1], samplestoplot[:, 2], color = colors[i + 1])
    end
    scatter!(ax, [mean(pcadf[:, colind1])], [mean(pcadf[:, colind2])], marker = :xcross,
             markersize = 20, color = :black)
    outliers = pcadf[abs.(pc1 .- meanpc1) .> outlierstd * stdpc1 .||
                     abs.(pc2 .- meanpc2) .> outlierstd * stdpc2, [1, colind1, colind2]]
    text!(ax, outliers[:, 1], position = [(row[2], row[3]) for row in eachrow(outliers)])
    println(outliers[:, 1])
    return (f, ax)
end

function graph3pcs(pcadf::AbstractDataFrame, outlierstd::Real = 2) # Needs GLMakie
    f = Figure()
    ax = Axis3(f[1,1])
    colors = [:blue, :yellow, :orange, :red]
    pc1 = pcadf[:, 3]
    pc2 = pcadf[:, 4]
    pc3 = pcadf[:, 5]
    meanpc1 = mean(pc1)
    meanpc2 = mean(pc2)
    meanpc3 = mean(pc3)
    stdpc1 = std(pc1)
    stdpc2 = std(pc2)
    stdpc3 = std(pc3)
    for i in 0:3
        samplestoplot = pcadf[abs.(pc1 .- meanpc1) .> i * stdpc1 .||
                              abs.(pc2 .- meanpc2) .> i * stdpc2 .||
                              abs.(pc3 .- meanpc3) .> i * stdpc3, 3:5]
        meshscatter!(ax, samplestoplot[:, 1], samplestoplot[:, 2], samplestoplot[:, 3],
                     markersize = 0.01, color = colors[i + 1])
    end
    meshscatter!(ax, [meanpc1], [meanpc2], [meanpc3],
             markersize = 0.01, color = :black)
    outliers = pcadf[abs.(pc1 .- meanpc1) .> outlierstd * stdpc1 .||
                     abs.(pc2 .- meanpc2) .> outlierstd * stdpc2 .||
                     abs.(pc3 .- meanpc3) .> outlierstd * stdpc3, [1; 3:5]]
    text!(ax, outliers[:, 1], position = [(row[2], row[3], row[4]) for row in eachrow(outliers)])
    println(outliers[:, 1])
    return f
end

pcadf = CSV.read(pcafile, DataFrame,
                 header = 0)
df = getoutliers(pcadf)
CSV.write(outfile * ".outlier", df, delim = ' ', writeheader = false)
for i in 1:5
    for j in (i+1):5
        f, ax = graphpcs(pcadf, i+2, j+2, 2)
        ax.xlabel = "PC$i"
        ax.ylabel = "PC$j"
        save(outfile * ".PC$(i)_PC$(j).png", f)
    end
end

