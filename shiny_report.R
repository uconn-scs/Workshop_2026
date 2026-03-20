library(shiny)
library(data.table)
library(DT)
library(msDiaLogue)
library(dplyr)
library(ggplot2)

options(shiny.sanitize.errors = FALSE)  # show real errors

ui <- fluidPage(
  
  titlePanel("msDiaLogue Shiny Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      
      fileInput(
        "specfile", "Upload Spectronaut CSV",
        accept = c(".csv", ".txt")
      ),
      
      checkboxInput("removeCont", "Remove Contaminants (CON__)", TRUE),
      
      actionButton("run", "Run Pipeline", class = "btn-primary"),
      
      hr(),
      
      uiOutput("ref_selector"),   # dynamic reference dropdown
      
      hr(),
      
      h4("Analysis Methods Used"),
      verbatimTextOutput("method_text")
      
    ),
    
    mainPanel(
      tabsetPanel(
        id = "tabs",
        tabPanel("Raw Histogram", plotOutput("hist_raw")),
        tabPanel("Normalization", plotOutput("norm_box")),
        tabPanel("Imputation", plotOutput("imp_dist")),
        tabPanel("Volcano", plotOutput("volcano")),
        tabPanel("PCA Score Plot", plotOutput("pca")),
        tabPanel("Summary Table", DTOutput("summary")),
        tabPanel("Fold Changes", DTOutput("fold_changes"))  # NEW TAB
      )
    )
  )
)



server <- function(input, output, session) {
  
  options(shiny.fullstacktrace = TRUE)
  
  #--------------------------------------------------------#
  # 0. Load raw file
  #--------------------------------------------------------#
  df_raw <- eventReactive(input$run, {
    req(input$specfile)
    
    df <- fread(input$specfile$datapath) |> as.data.frame()
    
    if (input$removeCont) {
      df <- df |> filter(!grepl("^CON__", PG.ProteinAccessions))
    }
    
    df
  })
  
  
  #--------------------------------------------------------#
  # 1. Preprocessing
  #--------------------------------------------------------#
  data_pre <- reactive({
    req(df_raw())
    
    preprocessing(
      dataSet = df_raw(),
      filterNaN = TRUE,
      filterUnique = 2,
      replaceBlank = TRUE,
      saveRm = FALSE
    )
  })
  
  
  #--------------------------------------------------------#
  # 2. Transformation (log2)
  #--------------------------------------------------------#
  data_trans <- reactive({
    req(data_pre())
    transform(data_pre(), method = "log", logFold = 2)
  })
  
  
  #--------------------------------------------------------#
  # 3. Normalization (sample-wise median)
  #--------------------------------------------------------#
  data_norm <- reactive({
    req(data_trans())
    normalize(data_trans(), 
              applyto = "sample", 
              normalizeType = "median", 
              plot = FALSE)
  })
  
  
  #--------------------------------------------------------#
  # 4. Filtering (minProp = 0.2)
  #--------------------------------------------------------#
  data_filt <- reactive({
    req(data_norm())
    filterNA(data_norm(), minProp = 0.2, by = "cond", saveRm = FALSE)
  })
  
  
  #--------------------------------------------------------#
  # 5. Imputation (seq-KNN)
  #--------------------------------------------------------#
  data_imp <- reactive({
    req(data_filt())
    impute.knn_seq(data_filt())
  })
  
  
  #--------------------------------------------------------#
  # 6. Dynamic Reference Selector
  #--------------------------------------------------------#
  output$ref_selector <- renderUI({
    req(data_imp())
    
    choices <- levels(data_imp()$R.Condition)
    
    selectInput(
      "ref_cond",
      "Select Reference Condition:",
      choices = choices,
      selected = choices[1]
    )
  })
  
  
  #--------------------------------------------------------#
  # 7. Moderated t-tests
  #--------------------------------------------------------#
  data_modt <- reactive({
    req(data_imp(), input$ref_cond)
    
    analyze.mod_t(
      data_imp(),
      ref = input$ref_cond,
      adjust.method = "BH"
    )
  })
  
  
  #--------------------------------------------------------#
  # 8. Methods description
  #--------------------------------------------------------#
  output$method_text <- renderText({
    paste(
      "• Log2 Transformation (base 2)\n",
      "• Median Sample-wise Normalization\n",
      "• seq-KNN Imputation (sequential nearest neighbors)\n",
      "• Moderated t-test (limma-style empirical Bayes shrinkage)\n",
      sep = ""
    )
  })
  
  
  #--------------------------------------------------------#
  # 9. Plots
  #--------------------------------------------------------#
  output$hist_raw <- renderPlot({
    req(df_raw())
    validate(need("PG.Quantity" %in% colnames(df_raw()), "Missing PG.Quantity column."))
    
    ggplot(df_raw(), aes(x = log2(PG.Quantity))) +
      geom_histogram() +
      ggtitle("Raw log2 Histogram")
  })
  
  
  output$norm_box <- renderPlot({
    req(data_norm())
    visualize.boxplot(data_norm())
  })
  
  
  output$imp_dist <- renderPlot({
    req(data_filt(), data_imp())
    visualize.dist(list("Before" = data_filt(), "After" = data_imp()))
  })
  
  
  output$volcano <- renderPlot({
    req(data_modt())
    
    firstContrast <- data_modt()[[1]]
    validate(need(!is.null(firstContrast), "No moderated t-test result available."))
    
    visualize.volcano(firstContrast)
  })
  
  
  output$pca <- renderPlot({
    req(data_imp())
    
    pca <- analyze.pca(data_imp())
    visualize.score(pca, ellipse = TRUE)
  })
  
  
  #--------------------------------------------------------#
  # 10. Summary Table
  #--------------------------------------------------------#
  output$summary <- renderDT({
    req(data_imp())
    
    summar <- summarize(data_imp(), saveSumm = FALSE) %>%
      mutate(across(where(is.numeric), ~ round(.x, 2)))
    
    datatable(summar, options = list(scrollX = TRUE, pageLength = 10))
  })
  
  
  #--------------------------------------------------------#
  # 11. NEW — Fold Change Table (all pairwise results)
  #--------------------------------------------------------#
  output$fold_changes <- renderDT({
    req(data_modt())
    
    modt_list <- data_modt()
    
    fc_table <- bind_rows(
      lapply(names(modt_list), function(comp) {
        
        df <- modt_list[[comp]]
        
        # Ensure expected rows exist
        validate(
          need("difference" %in% rownames(df),
               paste("No 'difference' row found for contrast:", comp))
        )
        validate(
          need("p-value" %in% rownames(df),
               paste("No 'p-value' row found for contrast:", comp))
        )
        
        # Extract values
        diffs <- df["difference", ]
        pvals <- df["p-value", ]
        
        # Build tidy format
        tibble(
          Comparison = comp,
          Protein = names(diffs),
          logFC = as.numeric(diffs),
          absFC = abs(as.numeric(diffs)),
          p_value = as.numeric(pvals)
        ) %>%
          mutate(
            adj_p_val = p.adjust(p_value, method = "BH")
          )
      })
    )
    
    # Sort
    fc_table <- fc_table %>% arrange(Comparison, desc(absFC))
    
    datatable(fc_table, options = list(scrollX = TRUE, pageLength = 15))
  })

}  
shinyApp(ui, server)