header <- dashboardHeader(title = "SharedPeptidesTMT")
sidebar <- dashboardSidebar(
    menuItem("Dataset view", tabName = "dataset"),
    menuItem("Cluster view: proteins", tabName = "cluster_prot"),
    menuItem("Cluster view: profiles", tabName = "cluster_profiles"),
    menuItem("Cluster view: what-if", tabName = "whatif"))
body <- dashboardBody(
    tabItems(
        tabItem(tabName = "dataset",
                fluidRow(DT::DTOutput("clusters_tbl"))),
        tabItem(tabName = "cluster_prot",
                fluidRow(DT::DTOutput("proteins_tbl"),
                         DT::DTOutput("proteins_summary_tbl"),
                         box(checkboxInput("use_unmodified", "Use unmodified peptides"),
                             checkboxInput("just_detectable", "Only use detectable peptides"),
                             checkboxInput("just_detected", "Only use detected peptides"),
                             plotOutput("prot_graph")))),
        tabItem(tabName = "cluster_profiles",
                fluidRow(
                    column(12,
                           box(width = 12,
                               selectInput("cluster_id",
                                           label = "Select a cluster",
                                           choices = unique(brd$Cluster),
                                           selected = 1),
                               checkboxInput("group_by_cond", "Group by Condition"),
                               checkboxInput("remove_control", "Remove control"),
                               checkboxInput("remove_outliers", "Remove outliers with PeCorA"),
                               checkboxInput("use_protein_ids", "Use protein IDs instead of gene IDs",
                                             value = TRUE),
                               sliderInput("pecora_threshold", "PeCorA p-value cutoff:",
                                           min = 0.01, max = 0.05, value = 0.01, step = 0.01),
                               selectInput("pecora_control", label = "PeCorA control group",
                                           choices = unique(brd$Condition))),
                           box(width = 12,
                               plotOutput("plot_cluster_full")))
                ),
                fluidRow(
                    column(12,
                           column(12,
                                  box("Summarization 1:",
                                      selectInput("summary_1",
                                                  label = "",
                                                  choices = c("All possible" = "all",
                                                              "Only unique" = "unique",
                                                              "Parsimony" = "parsimony",
                                                              "PeCorA" = "pecora",
                                                              "Quantifere" = "quantifere",
                                                              "PeCorA + Quantifere" = "pecora_quantifere")),
                                      plotOutput("summary_1_plot"),
                                      checkboxInput("sum1prof", "Keep feature profiles",
                                                    value = TRUE),
                                      width = 12)),
                           column(12,
                                  box("Summarization 2:",
                                      selectInput("summary_2",
                                                  label = "",
                                                  choices = c("All possible" = "all",
                                                              "Only unique" = "unique",
                                                              "Parsimony" = "parsimony",
                                                              "PeCorA" = "pecora",
                                                              "Quantifere" = "quantifere",
                                                              "PeCorA + Quantifere" = "pecora_quantifere")),
                                      plotOutput("summary_2_plot"),
                                      checkboxInput("sum2prof", "Keep feature profiles",
                                                    value = TRUE),
                                      width = 12)
                           )
                    )
                ),
                fluidRow(box("Quantifere clusters",
                             plotOutput("quantifere_groups"),
                             uiOutput("quantifere_selectors")
                ))

        ),
        tabItem(tabName = "whatif"))
)
dashboardPage(header, sidebar, body)

