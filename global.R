# Libraries ----
library(shiny)
library(MSstatsTMT)
library(shinydashboard)
library(SharedPeptidesExplorer)
library(data.table)
library(stringr)
library(PeCorA)
library(cleaver)
library(RcppGreedySetCover)
library(DT)
library(igraph)
# Data ----
brd = readRDS("www/brd.RDS")
brd_shared = brd[NumProteins > 1]
brd_prots_db = readRDS("www/just_brd.RDS")
detectability_data = readRDS("www/detectability_data.RDS")
# Functions ----
## PECORA
true_pecora_preprocessing = function(input, area_column_name = "Intensity",
                                     control_name = "PC9_DMSO") {
    t = as.data.frame(input)
    t["modpep_z"] = paste(t$Peptide.Modified.Sequence, "all",
                          sep = "_")
    t["condition_rep"] = paste(t$Condition, t$BioReplicate, sep = "_")
    if (suppressWarnings(length(which(is.na(as.numeric(t[, "Intensity"])) ==
                                      TRUE)) > 0)) {
        t <- t[-suppressWarnings(which(is.na(as.numeric(t[, area_column_name])) ==
                                           TRUE)), ]
    }
    t$ms1log2 <- log2(as.numeric(t[, area_column_name]))
    t$ms1scaled <- standardize::scale_by(ms1log2 ~ Condition * BioReplicate, t)
    ms1adj <- rep(0, nrow(t))
    i = 1
    t_cntrl <- t[as.character(t$Condition) == control_name, ]
    ms1scaled_cntrl <- t_cntrl[, "ms1scaled"]
    ms1scaled_full <- t[, "ms1scaled"]
    pb <- txtProgressBar(min = 0, max = length(unique(t$modpep_z)),
                         style = 3)
    for (x in unique(t$modpep_z)) {
        ms1adj[which(t$modpep_z == x)] <- suppressWarnings(ms1scaled_full[which(t$modpep_z ==
                                                                                    x)] - mean(ms1scaled_cntrl[which(t_cntrl$modpep_z ==
                                                                                                                         x)]))
        setTxtProgressBar(pb, i)
        i = i + 1
    }
    t$ms1adj <- ms1adj
    t
}
do_pecora_inference =  function(input, control_group) {
    pecora = unique(input[, .(Intensity = sum(Intensity)),
                          by = c("ProteinName", "PeptideSequence",
                                 "Condition", "BioReplicate")])
    setnames(pecora, "PeptideSequence", "Peptide.Modified.Sequence")
    processed = true_pecora_preprocessing(pecora, control_name = control_group)
    pecora_result = as.data.table(PeCorA(processed))
    pecora_result[, peptide := str_replace_all(peptide, "_all", "")]
    pecora_result[, is_outlier := TRUE]
    setnames(pecora_result,
             c("protein", "peptide"),
             c("ProteinName", "PeptideSequence"))
    function(cutoff) {
        outliers_1 = pecora_result[adj_pval <= cutoff]
        input = merge(input, outliers_1,
                      by = c("ProteinName", "PeptideSequence"),
                      all.x = TRUE, all.y = TRUE, allow.cartesian = TRUE)
        input[, is_outlier := ifelse(is.na(is_outlier), FALSE, is_outlier)]
        unique(input[!(is_outlier), .(ProteinName, PeptideSequence,
                                      Label = "PeCorA-adjusted")])
    }
}
## UNIQUE
do_unique_inference = function(input) {
    pp = unique(input[, .(ProteinName, PeptideSequence)])
    pp[, NumProteinsPerPeptide := uniqueN(ProteinName),
       by = "PeptideSequence"]
    function() {
        pp$Label = "Unique only"
        pp[NumProteinsPerPeptide == 1, .(ProteinName, PeptideSequence, Label)]
    }
}
## PARSIMONY
do_parsimony_inference = function(input) {
    pp = unique(input[, .(ProteinName, PeptideSequence)])
    pars = greedySetCover(pp)
    function() {
        pars$Label = "parsimony"
        pars
    }
}
## ALL
do_all_inference = function(input) {
    pp = unique(input[, .(ProteinName, PeptideSequence)])
    function() {
        pp$Label = "All possible"
        pp
    }
}
## QUANTIFERE
cluster_peptides = function(input, aggregate_bioreps = FALSE,
                            include_unique = FALSE,
                            use_control = FALSE, control_pattern = "DMSO") {
    pp = unique(input[, .(ProteinName, PeptideSequence)])
    if (include_unique) {
        pp_wide = dcast(unique(input[ ,
                                      .(ProteinName, PeptideSequence, Present = 1)]),
                        PeptideSequence ~ ProteinName, fill = 0)
    } else {
        pp_wide = dcast(unique(input[NumProteinsPerPeptide > 1,
                                     .(ProteinName, PeptideSequence, Present = 1)]),
                        PeptideSequence ~ ProteinName, fill = 0)

    }
    pp_wide2 = copy(pp_wide)
    pps_list = list()
    not_end = TRUE
    while (not_end) {
        peptide = unique(pp_wide2$PeptideSequence)[1]
        assignments = unlist(pp_wide[PeptideSequence == peptide, -1])
        assignments = assignments[assignments == 1]
        proteins = names(assignments)
        pps = unique(pp[ProteinName %in% proteins, PeptideSequence])
        pp_wide2 = pp_wide2[!(PeptideSequence %in% pps)]
        pps_list = c(pps_list, list(pps))
        if (nrow(pp_wide2) == 0) {
            not_end = FALSE
        }
    }

    pps_df = rbindlist(lapply(seq_along(pps_list),
                              function(x) data.table(PeptideGroup = x,
                                                     PeptideSequence = pps_list[[x]])))
    all_pp = merge(pp, pps_df, by = "PeptideSequence")
    with_group = merge(input, all_pp, by = c("ProteinName", "PeptideSequence"),
                       allow.cartesian = TRUE)
    with_group_split = split(with_group, with_group$PeptideGroup)
    clustering = lapply(with_group_split, function(x) {
        just_features = unique(x[,
                                 .(PeptideSequence, Charge, PSM, Condition, BioReplicate,
                                   Run, TechRepMixture, Mixture, Channel, Intensity)])
        if (!use_control) {
            just_features = just_features[!grepl(control_pattern, BioReplicate, fixed = TRUE)]
        }
        if (aggregate_bioreps) {
            just_features = just_features[,
                                          .(Intensity = mean(Intensity,
                                                             na.rm = TRUE)),
                                          setdiff(colnames(just_features),
                                                  c("BioReplicate", "Intensity",
                                                    "Channel"))]
            just_features[, ProfileColumn := Condition]
        } else {
            just_features[, ProfileColumn := BioReplicate]
        }
        n_items = uniqueN(just_features$PeptideSequence)
        just_features[, log2Intensity := log(Intensity, 2)]
        just_features = just_features[, .(PSM, ProfileColumn, log2Intensity)]
        just_features[, log2Intensity := scale(log2Intensity),
                      by = "PSM"]
        features_wide = dcast(just_features, ProfileColumn ~ PSM)
        features_wide = features_wide[, -1, with = FALSE]
        dists = proxy::dist(t(features_wide), method = function(x, y) 1 - cor(x, y))
        clusters = hclust(dists, method = "average")
        clusters
    })
    clustering
}
do_quantifere_inference = function(input, clustering) {
    pp = unique(input[, .(ProteinName, PeptideSequence)])
    primary_proteins = do_parsimony_inference(pp)()
    primary_protein_ids = unique(primary_proteins$ProteinName)

    function(num_clusters) {
        peptide_groups = lapply(1:length(clustering), function(i) {
            groups = cutree(clustering[[i]], k = num_clusters[i])
            peptide_groups = data.table(PSM = names(groups),
                                        PeptideCluster = groups)
            peptide_groups$PeptideCluster = paste(i,
                                                  peptide_groups$PeptideCluster,
                                                  sep = "_")
            peptide_groups

        })
        peptide_groups = rbindlist(peptide_groups)
        full_cluster_info = merge(peptide_groups, unique(input[, .(ProteinName, PSM)]),
                                  by = "PSM", all.x = TRUE,
                                  all.y = TRUE, allow.cartesian = TRUE)
        full_cluster_info$IsPrimary = full_cluster_info$ProteinName %in% primary_protein_ids
        by_group = split(full_cluster_info, full_cluster_info$PeptideCluster)
        secondary_proteins = lapply(by_group, function(group_df) {
            if (all(!group_df$IsPrimary)) {
                unique(group_df$ProteinName)
            } else {
                character(0)
            }
        })
        secondary_proteins = unlist(secondary_proteins, use.names = FALSE)
        pp[ProteinName %in% c(primary_protein_ids, secondary_proteins),
           .(ProteinName, PeptideSequence, Label = "Quantifere")]
    }
}
## PECORA+QUANTIFERE
# do_pecora_quantifere = function(input, pecora_inferred) {
#
# }

# App ----
# Cluster view
# unique(brd[, .(Cluster, NumProteins, NumPeptides, AnyHasUnique,
#                EachHasUnique)])

# get_protein_summary = function(input, inferred_proteins, ...) {
#     inferred = do.call("inferred_proteins", args = list(...))
#     input = merge(input, inferred,
#                   by = c("ProteinName", "PeptideSequence"))
#     by_protein = split(input, input$ProteinName)
#     summarized = lapply(by_protein, MSstatsTMT::proteinSummarization)
#     rbindlist(summarized)
# }

update_profiles = function(profiles, inferred) {
    inferred$IsInferred = TRUE
    input = merge(profiles, inferred,
                  by = c("ProteinName", "PeptideSequence"),
                  all.x = TRUE)
    input$IsInferred = ifelse(is.na(input$IsInferred), FALSE, input$IsInferred)
    input$Label = unique(inferred$Label, na.rm = TRUE)
    input
}

get_protein_summary = function(updated_profiles, inferred, log_trans = FALSE) {
    input = updated_profiles[(IsInferred)]
    if (!log_trans) {
        input$Intensity = 2^input$Intensity
    }
    by_protein = split(input, input$ProteinName)
    summarized = lapply(by_protein,
                        function(x) MSstatsTMT::proteinSummarization(unique(x),
                                                                     global_norm = FALSE))
    summarized = rbindlist(summarized)
    summarized = merge(summarized, unique(updated_profiles[, .(GeneName, ProteinName)]),
                       by.x = "Protein", by.y = "ProteinName")
    summarized$Label = unique(inferred$Label)
    summarized
}

get_summary_comparison = function(summary_1, summary_2) {
    comp = rbind(summary_1, summary_2)
    comp
}

plot_summary = function(profiles, summary, profile_column = "BioReplicate",
                        plot_profiles = TRUE) {
    plot = ggplot(profiles,
                  aes_string(x = profile_column, y = "Intensity",
                             group = "PSM", linetype = "IsInferred")) +
        geom_line(aes_string(x = profile_column, y = "Abundance", group = "Label",
                             color = "Label"),
                  data = summary,
                  inherit.aes = FALSE, size = 2) +
        facet_wrap(~Group) +
        theme_bw()
    if (plot_profiles) {
        plot = plot + geom_line(color = "grey", size = 1.2)
    }
    plot
}


do_protein_inference = function(input, method, clusters = NULL,
                                control = NULL, pecora_cutoff = NULL,
                                quantifere_clusters = NULL) {
    pp = unique(input[, .(ProteinName, PeptideSequence)])
    if (method == "parsimony") {
        result = do_parsimony_inference(input)()
    } else if (method == "unique") {
        result = do_unique_inference(input)()
    } else if (method == "all") {
        result = do_all_inference(input)()
    } else if (method == "pecora") {
        result = do_pecora_inference(input, control)(pecora_cutoff)
    } else if (method == "quantifere") {
        result = do_quantifere_inference(input, clusters)(quantifere_clusters)
    }
    result
}

plot_summary = function(profiles, summary,
                        profile_column = "BioReplicate",
                        plot_profiles = TRUE) {
    profiles$IsInferred = factor(as.character(ifelse(profiles$IsInferred, "yes", "no")),
                                 levels = c("yes", "no"), ordered = TRUE)
    setnames(summary, "Protein", "ProteinName")
    plot = ggplot(profiles,
                  aes_string(x = profile_column, y = "Intensity",
                             group = "PSM", linetype = "IsInferred")) +
        theme_bw()
    if (plot_profiles) {
        plot = plot + geom_line(color = "grey", size = 1.5)
    }
    plot + geom_line(aes_string(x = profile_column, y = "Abundance", group = "Label",
                                color = "Label"),
                     data = summary,
                     inherit.aes = FALSE, size = 2) +
        facet_wrap(~Group) +
        theme(axis.text.x = element_text(angle = 270))
}


plot_quantifere_groups = function(input, profile_column = "BioReplicate") {
    ggplot(input, aes_string(x = profile_column, y = "Intensity",
                             group = "PSM", color = "PeptideSequence")) +
        geom_point() +
        geom_line() +
        facet_wrap(~Group) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 270)) +
        theme(legend.position = "bottom")
}


get_protein_graph = function(input, cluster_id, use_unmodified, detectable_only,
                                 detected_only) {
    if (use_unmodified) {
        cols = c("ProteinName", "UnmodifiedSequence")
    } else {
        cols = c("ProteinName", "PeptideSequence")
    }
    if (detectable_only) {
        input = input[(IsDetectable) | (IsDetected)]
    }
    if (detected_only) {
        input = input[(IsDetected)]
    }
    to_graph = unique(input[Cluster == cluster_id, cols, with = FALSE])
    setnames(to_graph, cols, c("ProteinName", "PeptideSequence"))
    to_graph$Present = 1
    wide = dcast(to_graph, ProteinName ~ PeptideSequence, fill = 0)
    wide_m = as.matrix(wide[, -1, with = FALSE])
    row.names(wide_m) = wide[, ProteinName]
    pp_g = wide_m %*% t(wide_m)
    diag(pp_g) = 0
    counts = to_graph[, .(NumPeptides = uniqueN(PeptideSequence)),
                      by = "ProteinName"]
    g = graph_from_adjacency_matrix(pp_g,
                                    weighted = TRUE, mode = "undirected")
    plot(g, edge.label = E(g)$weight, vertex.label = paste(names(V(g)), counts$NumPeptides,
                                                           sep = " - "))
}
