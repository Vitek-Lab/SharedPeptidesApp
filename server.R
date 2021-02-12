server <- function(input, output) {

    grouping_column = reactive({
        if (input$use_protein_ids) {
            "ProteinName"
        } else {
            "GeneName"
        }
    })

    current_cluster = reactive({
        current = brd[Cluster == input$cluster_id]
        if (input$remove_control) {
            current = current[!grepl("DMSO", BioReplicate)]
        }
        if (input$group_by_cond) {
            current = current[, .(Intensity = mean(Intensity, na.rm = TRUE)),
                              by = setdiff(colnames(current), c("BioReplicate",
                                                                "Intensity",
                                                                "Channel"))]
        }
        current$Intensity = log(current$Intensity, 2)
        if (input$remove_outliers) {
            inferred = do_pecora_inference(current, input$pecora_control)
            inferred_df = inferred(input$pecora_threshold)
            current = merge(current, inferred_df[, .(ProteinName, PeptideSequence)],
                            by = c("ProteinName", "PeptideSequence"))
        }
        current$Group = current[[grouping_column()]]
        current
    })

    output$plot_cluster_full = renderPlot({
        if (input$group_by_cond) {
            x_axis = "Condition"
        } else {
            x_axis = "BioReplicate"
        }
        plot = ggplot(current_cluster(), aes_string(x = x_axis, y = "Intensity",
                                             group = "PSM")) +
            # geom_line() +
            facet_wrap(~Group) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 270))
        if (!input$use_protein_ids) {
            plot = plot +
                geom_line(aes(color = ProteinName), size = 1.2) +
                theme(legend.position = "bottom")
        } else {
            plot = plot +
                geom_line(size = 1.2)
        }
        plot
    })

    clusters_quantifere = reactive({
        if (input$summary_2 == "quantifere" | input$summary_1 == "quantifere") {
            clusters = cluster_peptides(current_cluster(), input$group_by_cond, TRUE)
        } else {
            clusters = NULL
        }
        clusters
    })

    inference_1 = reactive({
        clusters = clusters_quantifere()
        inference = do_protein_inference(current_cluster(), input$summary_1,
                                         clusters,
                                         control = input$pecora_control,
                                         pecora_cutoff = input$pecora_threshold,
                                         quantifere_clusters = quantifere_nums())
        inference
    })
    inference_2 = reactive({
        clusters = clusters_quantifere()
        inference = do_protein_inference(current_cluster(), input$summary_2,
                                         clusters,
                                         control = input$pecora_control,
                                         pecora_cutoff = input$pecora_threshold,
                                         quantifere_clusters = quantifere_nums())
        inference
    })


    output$summary_1_plot = renderPlot({
        updated_profiles = update_profiles(current_cluster(), inference_1())
        summary = get_protein_summary(updated_profiles, inference_1())
        plot_summary(updated_profiles, summary, plot_profiles = input$sum1prof)
    })
    output$summary_2_plot = renderPlot({
        updated_profiles = update_profiles(current_cluster(), inference_2())
        summary = get_protein_summary(updated_profiles, inference_2())
        plot_summary(updated_profiles, summary, plot_profiles = input$sum2prof)
    })

    quantifere_nums = reactive({
        if (input$summary_1 == "quantifere" | input$summary_2 == "quantifere") {
            nums = vector("integer", length(clusters_quantifere()))
            for (i in 1:length(clusters_quantifere())) {
                nums[i] = input[[paste0("num_quantifere", i)]]
            }
            nums
        } else {
            NULL
        }
    })

    output$quantifere_groups = renderPlot({
        clusters = clusters_quantifere()

        feature_groups = rbindlist(lapply(seq_along(clusters), function(i) {
            data.table(PeptideGroup = i,
                       PSM = clusters[[i]]$labels)
        }))

        for_group_plot = merge(current_cluster(), feature_groups, by = "PSM",
                               all.x = TRUE, all.y = TRUE,
                               allow.cartesian = TRUE)
        plot_quantifere_groups(for_group_plot,
                               ifelse(input$group_by_cond,
                                      "Condition", "BioReplicate"))
    })

    output$quantifere_selectors = renderUI({
        n = length(clusters_quantifere())
        output_list <- lapply(1:n, function(i) {
            name <- paste("select_quantifere", i, sep="")
            uiOutput(name)
        })

        # Convert the list to a tagList - this is necessary for the list of items
        # to display properly.
        do.call(tagList, output_list)
    })

    observe({
        for (i in 1:length(clusters_quantifere())) {
            local({
                my_i <- i
                name <- paste("select_quantifere", my_i, sep = "")

                output[[name]] <- renderUI({
                    selectInput(paste0("num_quantifere", my_i),
                                "Number of peptide clusters in the group:",
                                choices = 1:length(clusters_quantifere()[[my_i]]$labels))
                })
            })
        }
    })

    output$clusters_tbl = DT::renderDT({
        unique(brd[, .(Cluster, NumProteins, NumPeptides, AnyHasUnique,
                       EachHasUnique)])
    }, filter = "top")


    output$proteins_tbl = DT::renderDT({
        unique(detectability_data[Cluster == input$cluster_id,
                                  .(GeneName, ProteinName, UnmodifiedSequence, PeptideSequence,
                                    IsDetected, IsTryptic, IsModified,
                                    IsDetectable, DetectionProbability, NumVariants,
                                    AlgoFailed)])
    }, filter = "top")


    output$proteins_summary_tbl = DT::renderDT({
        result = detectability_data[Cluster == input$cluster_id]
        result[, IsUnique := uniqueN(ProteinName) == 1, by = "UnmodifiedSequence"]
        result[, IsUniqueObserved := (uniqueN(ProteinName) == 1), by = "PeptideSequence"]
        result = result[, .(NumDetected = sum(IsDetected, na.rm = TRUE),
                            HasUniqueObserved = any(IsUniqueObserved),
                            NumNotDetected = sum(!IsDetected & IsDetectable, na.rm = TRUE),
                            TotalDetectable = sum(IsDetectable, na.rm = TRUE),
                            NumModified = sum(IsModified, na.rm = TRUE),
                            UniqueDetectable = sum(IsUnique & IsDetectable, na.rm = TRUE)),
                        by = "ProteinName"]
        result
    })

    output$prot_graph = renderPlot({
        get_protein_graph(detectability_data, input$cluster_id,
                          input$use_unmodified, input$just_detectable,
                          input$just_detected)
    })
}
