library(shiny) #version 1.7.4
library(shinyFiles) #version 0.9.3
library(shinycssloaders) #version 1.0.0
library(shinyWidgets) #version 0.8.0
library(data.table) #version 1.14.2
library(ggplot2) #version 3.3.5
library(dplyr) #version 1.0.7
library(plyr) #version 1.8.6
library(tidyr) #version 1.1.4
library(naniar) #version 1.0.0
library(gridExtra) #version 2.3
library(tidyverse) #version
library(grid)#version 4.1.2
library(DT)#version 0.26
library(ggrepel) #version 0.9.3

ui <- fluidPage(
  
  titlePanel("Sex Identification"),
  
  sidebarPanel(
    wellPanel(
      tags$h4("Upload your data:"),
      fileInput("table1", label = h5("Peptide Intensities Data")),
      fileInput("table2", label = h5("Experimental Design")),
      tags$h5("Filtering parameters:"),
      numericInput("dotProduct", "Library Dot Product", value = 0.9, min = 0, max = 1, step = 0.05),
      numericInput("isotopeDotProduct", "Isotope Dot Product", value = 0.9, min = 0, max = 1, step = 0.05),
      numericInput("peakFoundRatio", "Peptide Peak Found Ratio", value = 0.9, min = 0, max = 1, step = 0.05),
      actionButton("stdcurvesButton", "Plot STD curves")
    ),
    wellPanel(
      tags$h4("Model Assessment"),
      sliderInput("outlier", "max % of data removed", value = 0, min = 0, max = 0.5, step=0.05), 
      radioButtons("model",
                   h5("Model options"), 
                   choices = list("Experimental model" = 1, 
                                  "Pre-defined model" = 2),
                   selected = 1),
      actionButton("modelButton", "Plot model")
    ),
    wellPanel(
      tags$h4("Sex Identification"),
      actionButton("sexidButton", "Plot table"),
      downloadButton('download',"Download the data"),
    ),
    wellPanel(
      tags$h4("Sample information"),
      fileInput("skyline", label = h5("Raw Intensities Data ")),
      uiOutput("samples"),
      actionButton("tracesButton", "Plot XIC")
    )
  ),
  
  mainPanel(
    tabsetPanel(type = "tabs",
                tabPanel("Standards",mainPanel(h3(""), plotOutput("plot1", height="1000px", width = "1400px") %>% withSpinner(color="#0dc5c1"))
                ),
                tabPanel("Summary",
                         mainPanel(h3("Model"), plotOutput("plot2", width = "1400px") %>% withSpinner(color="#0dc5c1")),
                         mainPanel(h3("Results"), DT::dataTableOutput("result", width = "1400px") %>% withSpinner(color="#0dc5c1")),
                ),
                tabPanel("Signal per sample",
                         mainPanel(h3("MS2 traces"), plotOutput("plot3", height = "1000px", width = "1400px") %>% withSpinner(color="#0dc5c1")),
                )
    )
  )
)




server <- (function(input, output, session) {
  options(shiny.maxRequestSize=90*1024^2)
  

################################################################################################################################################
  #0. Load experimental design table - (fn = Experimental design)
  load_ExpDesign <- function(fn){
    ExpDesign <- fread(fn)
    Names <- c("File Name", "Experiment")
    Col_Names <- colnames(ExpDesign)
    ifelse(identical(Names,Col_Names) == FALSE, stop("Error"), "")
    return(ExpDesign)
  }
  
  #1. Load peptide intensities table - (fn = peptide intensities)
  load_table <- function (fn){
    Table <- fread(fn)
    #Check that the table is correct 
    Names <- c("Peptide", "Protein", "Replicate", "Peptide Peak Found Ratio",
               "Normalized Area", "Library Dot Product", "Isotope Dot Product",
               "Peptide Modified Sequence", "File Name", "Sample Type", "Analyte Concentration", "Concentration Multiplier")
    Col_Names <- colnames(Table)
    ifelse(identical(Names,Col_Names) == FALSE, stop("Error"), "")
    print(Table)
    return(Table)
  }
  
  #2. Filter peptide intensity table based on idopt, dopt and ppfr - (dt = Table)
  Filter <- function(dt, dotProd, isoDotProd, peakRatio) {
    dt$`Analyte Concentration` <-dt$`Analyte Concentration`*dt$`Concentration Multiplier`
    dt <- dt[dt$`Library Dot Product` >= dotProd,]
    dt <- dt[dt$`Isotope Dot Product` >= isoDotProd,]
    dt <- dt[dt$`Peptide Peak Found Ratio` >= peakRatio,]
    return(dt)
  }
  
  #3. Load and filter tables -  (fn1 = Experimental design, fn2 = Peptide intensities)
  Load_Filter <- function (fn1, fn2, dotProd, isoDotProd, peakRatio) {
    ExpDesign <- load_ExpDesign(fn2)
    Table <- load_table(fn1)
    Table <- inner_join(Table, ExpDesign , by="File Name", relationship = "many-to-many")
    TableF <- Filter(Table, dotProd, isoDotProd, peakRatio) #Filter table based on PPFR, IDOP, DOP
    return(TableF)
  }
  
  #4. STD table - Subset the standards from the table - (dt = TableF)
  #Input: Filtered table. Output: Table with STD only
  STD_Table <- function(dt){
    STD <- dt [dt$`Sample Type` == "Standard", ]
    STD <- STD[!is.na(STD$`Normalized Area`), ]
    colnames(STD) <- c(
      "Peptide",
      "Protein",
      "Replicate",
      "Peptide Peak Found Ratio",
      "Normalized Area",
      "Library Dot Product",
      "Isotope Dot Product",
      "Sequence",
      "File Name",
      "Sample Type",
      "Analyte Concentration",
      "Concentration Multiplier",
      "Experiment")
    return(STD)
  }
  
  #5. Regression - Perform linear regression and store the relevant coefficients - (dt = Table_i) 
  #Input: Table with standards for 1 experiment. Output: list with regression coefficients
  Regression <- function(dt) {
    Reg <-lm(dt$`Analyte Concentration` ~ dt$`Normalized Area`, data = dt)
    Slope <- round(coef(Reg)[2], 15)
    Intercept <- round(coef(Reg)[1], 3)
    R2 <- round(as.numeric(summary(Reg)[8]), 3)
    SDy <- as.numeric(summary(Reg)[["coefficients"]][1,2])
    LOD <- (3.3 * SDy) / Slope
    LOQ <- (10 * SDy) / Slope
    c(Slope, Intercept, R2, LOD, LOQ)
  }
  
  #6. Perform regression on the Std - (dt = STD)
  #Input: Table with STD . Output: Table with regression values
  STD_Regression <- function(dt){
    Reg <- data.table(Sequence = character(), Slope = numeric(), Intercept = numeric(), R2 = numeric(), LOD = numeric(), LOQ = numeric(), Experiment = character())
    for (i in 1:length(unique(unlist(dt[,Experiment])))) {
      Table_i <- dt[Experiment == unique(unlist(dt[,Experiment]))[i], ]
      regressions_data <- as.data.table(plyr::ddply(Table_i, "Sequence", Regression))
      regressions_data <- regressions_data[, Experiment := unique(unlist(dt[,Experiment]))[i]]
      colnames(regressions_data) <-
        c ("Sequence", "Slope", "Intercept", "R2", "LOD", "LOQ", "Experiment")
      Reg <- rbind(Reg, regressions_data)
    }
    return(Reg)
  }
  
  #7. Regression 2 - Generate the regression for the male model with the outlier control - (dt = Males, a = max % of values to be removed)
  #Inputs: Table with males only, number between 0 and 1 representing the max % of data that can be removed to build the model. Output: list with regression coefficients and filtered tables with the outliers removed. 
  Regression2 <- function(dt, a) {
    outliers_removed <- 0
    while (outliers_removed < a * nrow(dt)) {
      # Perform  linear regression
      Reg <- lm(log(dt$`AMELY`, 2) ~ log(dt$`AMELX`, 2), data = dt)
      # Calculate residuals from the regression
      residuals <- residuals(Reg)
      # Calculate standard deviation of residuals
      sd_residuals <- sd(residuals)
      # Identify the index of the highest residual
      max_residual_index <- which.max(abs(residuals))
      # Check if the highest residual is outlier
      if (abs(residuals[max_residual_index]) > 2 * sd_residuals) {
        dt <- dt[-max_residual_index, ]
        outliers_removed <- outliers_removed + 1
      } else {
        break  # Stop the loop if all points are within 2 SD
      }
    }
    # Perform final robust linear regression on the filtered data
    Reg <- lm(log(dt$`AMELY`, 2) ~ log(dt$`AMELX`, 2), data = dt)
    Slope <- round(Reg$coefficients[2], 3)
    Intercept <- round(Reg$coefficients[1], 3)
    R2 <- round(as.numeric(summary(Reg)[8]), 3)
    return(list(ModelInfo = list(Slope = Slope, Intercept = Intercept, R2 = R2), FilteredData = dt))
  }
  
  #8. Get Male identification - Filter when intensity < LOD and identify samples with AMELY as males - (dt1 =  tableF, dt2 = Reg)
  #Input: Filtered intensity table, Table with regression coefficients from the standards. Output: Table annotated with samples identified as males. 
  Male_ID <- function(dt1, dt2){
    Sample <- dt1[dt1$`Sample Type` == "Unknown",]
    NewTable <- data.table()
    
    # Loop through each experiment safely
    for (exp in unique(Sample$Experiment)) {
      data <- Sample[Sample$Experiment == exp, ]
      Stat <- dt2[dt2$Experiment == exp, ] # This will be an empty table if no standards for this exp
      
      datawide <- LOQ_filter(data, Stat)
      
      if(nrow(datawide) > 0) {
        datawide[, Experiment := exp]
        NewTable <- rbind(NewTable, datawide, fill = TRUE)
      }
    }
    
    # Handle case where NewTable is completely empty
    if(nrow(NewTable) == 0) {
      return(data.table(Replicate = character(), AMELX = numeric(), AMELY = numeric(), Sex = character()))
    }
    
    # Add sex identification
    # Ensure AMELX and AMELY columns exist to avoid errors
    if (!"AMELY" %in% names(NewTable)) NewTable[, AMELY := 0]
    if (!"AMELX" %in% names(NewTable)) NewTable[, AMELX := 0]
    
    NewTable[, Sex := fifelse(AMELY == 0, "Unknown", 
                              fifelse(AMELX == 0, "NonConclusive", "Male"))]
    return(NewTable)
  }
  
  #8.2. Filter out values that are < LOQ per target for a given experiment 
  LOQ_filter <- function(data, Stat){
    data <- data[, c("Replicate", "Peptide Modified Sequence", "Normalized Area")]
    data <- na.omit(data)
    
    if (nrow(data) == 0) {
      return(data.table(Replicate=character(), AMELX=numeric(), AMELY=numeric()))
    }
    
    datawide <- as.data.table(spread(data, "Peptide Modified Sequence", "Normalized Area"))
    
    # Create a named vector for LOQs for easy and safe lookup.
    # If Stat is empty or has no LOQ column, this will be an empty named vector.
    loq_lookup <- if(nrow(Stat) > 0 && "LOQ" %in% names(Stat)) setNames(Stat$LOQ, Stat$Sequence) else setNames(numeric(), character())
    
    # Helper function to get LOQ for a peptide, defaulting to 0 if not found.
    # A LOQ of 0 means the filter will never be applied if the standard is missing.
    get_loq <- function(peptide_name) {
      val <- loq_lookup[peptide_name]
      if (is.na(val) || is.null(val)) 0 else val
    }
    
    # Apply filter if the column exists
    if ("SIRPPYPSY" %in% names(datawide)) datawide$SIRPPYPSY[datawide$SIRPPYPSY < get_loq("SIRPPYPSY")] <- 0
    if ("SIRPPYPSYG" %in% names(datawide)) datawide$SIRPPYPSYG[datawide$SIRPPYPSYG < get_loq("SIRPPYPSYG")] <- 0
    if ("SM[+16]IRPPY" %in% names(datawide)) datawide$`SM[+16]IRPPY`[datawide$`SM[+16]IRPPY` < get_loq("SM[+16]IRPPY")] <- 0
    if ("SM[+16]IRPPYS" %in% names(datawide)) datawide$`SM[+16]IRPPYS`[datawide$`SM[+16]IRPPYS` < get_loq("SM[+16]IRPPYS")] <- 0
    if ("SMIRPPY" %in% names(datawide)) datawide$SMIRPPY[datawide$SMIRPPY < get_loq("SMIRPPY")] <- 0
    
    # Turns NA into 0 for all columns
    for (j in seq_along(datawide)) {
      set(datawide, which(is.na(datawide[[j]])), j, 0)
    }
    
    # Summing up peptides. Checks if columns exist before summing.
    amelx_peps <- c("SIRPPYPSY", "SIRPPYPSYG")
    amely_peps <- c("SM[+16]IRPPY", "SM[+16]IRPPYS", "SMIRPPY")
    
    datawide[, AMELX := rowSums(.SD, na.rm = TRUE), .SDcols = intersect(amelx_peps, names(datawide))]
    datawide[, AMELY := rowSums(.SD, na.rm = TRUE), .SDcols = intersect(amely_peps, names(datawide))]
    
    return(datawide)
  }
  
  #9.1. Male model predefined - Calculate ThAMELY (predicted AMELY intensity) of samples where no AMELY could be quantified with confidence intervals based on pre-defined model.- (dt = ID_Male)
  #Input: output from Male_ID function. Output: Table with predicted AMELY intensity for potential females with 95% confidence interval.
  Male_Predef_Model <- function(dt) {
    #Generate table with model coefficients
    coefs <- c("a", "b", "c", "r2")
    line <- c("0.7401", "8.0835", "NA", "0.9118")
    up <- c("0.0218", "-0.418", "23.6", "1")
    down <- c("-0.0218", "1.9", "-7.4", "1")
    Model <- data.table(coefs, line, up, down)
    #For AMELX in females calculate ThAMELY and AMELY_Low and AMELY_high
    #Subset the potential females
    Unknown <- dt[dt$Sex == "Unknown",]
    #Calculate the theoretical AMELY intensity based on the linear model
    Unknown [, ThAMELY := (as.numeric(Model[1, 2]) * log(AMELX, 2) + as.numeric(Model[2, 2]))]
    #Calculate the upper and lower AMELY values based on the confidence interval of the model
    Unknown [, AMELY_low := (as.numeric(Model[1, 4]) * log(AMELX,2) * log(AMELX,2)) + (as.numeric(Model[2, 4]) * log(AMELX, 2)) + as.numeric(Model[3, 4])]
    Unknown [, AMELY_high := (as.numeric(Model[1, 3]) * log(AMELX, 2) * log(AMELX, 2)) + (as.numeric(Model[2, 3]) * log(AMELX, 2)) + as.numeric(Model[3, 3])]
    return(Unknown)
  }
  
  #9.2. Generate the male model - calculate the predicted AMELY intensity with 95% confidence interval for potential females - (dt = ID_Male, a = max % of values to be removed)
  #Input: output from Male_ID function, number between 0 and 1 representing the max % of data that can be removed to build the model. Output: Table with predicted AMELY intensity for potential females with 95% confidence interval.
  Male_model <- function(dt, a) {
    #Subset the males to generate model
    Males <- dt[dt$Sex == "Male", ]
    #generate linear model of AMELX ~ AMELY intensity
    regressions_data2 <- Regression2(Males, a)
    #Get the confidence intervals of the male model and model the min and max at 95% confidence (ie. 2*error)
    FilteredMale <- as.data.table(regressions_data2$FilteredData)
    smooth_values <-
      predict(lm(log(AMELY, 2) ~ log(AMELX,2), FilteredMale),
              data = log(FilteredMale$AMELX,2),
              se.fit = TRUE)
    smooth <- as.data.table(smooth_values)
    ModelInfo <-regressions_data2$ModelInfo
    FilteredMale [, Line := ModelInfo$Slope * log(AMELX,2) + ModelInfo$Intercept]
    FilteredMale [, y_sepos := Line + 2 * smooth$se.fit]
    FilteredMale [, y_seneg := Line - 2 * smooth$se.fit]
    #Model with a polynomial function of degree 2 the upper and lower limits of the confidence interval and get the equation coefficients (y=ax2+bx+c)
    Model_up <-
      lm(FilteredMale$y_sepos ~ poly(log(FilteredMale$AMELX,2), 2, raw = TRUE), data = FilteredMale)
    a_up <- coef(Model_up)[3]
    b_up <- coef(Model_up)[2]
    c_up <- coef(Model_up)[1]
    Model_down <-
      lm(FilteredMale$y_seneg ~ poly(log(FilteredMale$AMELX,2), 2, raw = TRUE), data = FilteredMale)
    a_down <- coef(Model_down)[3]
    b_down <- coef(Model_down)[2]
    c_down <- coef(Model_down)[1]
    #Subset the potential females
    Unknown <- dt[dt$Sex == "Unknown", ]
    #Calculate the theoretical AMELY intensity based on the linear model
    Unknown [, ThAMELY := ModelInfo$Slope * log(AMELX,2) + ModelInfo$Intercept]
    #Calculate the upper and lower AMELY values based on the confidence interval of the model
    Unknown [, AMELY_low := (a_down * log(AMELX, 2) * log(AMELX, 2)) + (b_down * log(AMELX,2)) + c_down]
    Unknown [, AMELY_high := (a_up * log(AMELX,2) * log(AMELX,2)) + (b_up * log(AMELX, 2)) + c_up]
    return(Unknown)}
  
  #10. Get female identification - Confirm the female identification when ThAMELY fit predefined criteria. - (dt1 = AMELYTh , dt2 = TableF)
  #Input: Table with predicted AMELY intensities for potential females, Peptide intensities filtered table. Output: Table with female identification.
  Female_ID <- function(dt1, dt2){
    STD <- STD_Table(dt2)
    Reg <- STD_Regression(STD)
    ID <- data.table() # Initialize an empty data.table
    
    # A small helper function to safely get max value, returns -Inf if no peptides are found
    safe_max <- function(values) {
      if (length(values) == 0 || all(is.na(values))) return(-Inf)
      max(values, na.rm = TRUE)
    }
    
    for (exp in unique(dt2$Experiment)) {
      Reg2 <- Reg[Reg$Experiment == exp, ]
      Table <- dt1[dt1$Experiment == exp, ]
      
      # If there's no regression data for this experiment, mark all as NonConclusive
      if (nrow(Reg2) == 0) {
        Table[, Females := "NonConclusive"]
        ID <- rbind(ID, Table, fill = TRUE) # Use fill=TRUE to handle potential column mismatches
        next # Skip to the next iteration of the loop
      }
      
      # Robustly find LOD and LOQ by peptide name
      LOD_AMELY_val <- safe_max(Reg2[Sequence %in% c("SM[+16]IRPPY", "SM[+16]IRPPYS", "SMIRPPY"), LOD])
      LOQ_AMELY_val <- safe_max(Reg2[Sequence %in% c("SM[+16]IRPPY", "SM[+16]IRPPYS", "SMIRPPY"), LOQ])
      
      # Apply the logic
      Table[, Females := fifelse(ThAMELY > log(LOQ_AMELY_val, 2) & AMELY_low > log(LOD_AMELY_val, 2), "Female", "NonConclusive")]
      ID <- rbind(ID, Table, fill = TRUE)
    }
    # Ensure the final columns match the expected output, filling missing ones with NA
    # This prevents errors if 'Females' column isn't created in some loops
    if (!"Females" %in% names(ID)) {
      ID[, Females := NA_character_]
    }
    return(ID)
  }
  #11. Combine identifications in summary table - (dt1 = Female, dt2 = ID_Male)
  #Input: Table with female identification, Table with male identification. Output: Combined table with both male and female identification.
  Summary_table <- function (dt1, dt2) {
    Males2 <- dt2 [dt2$Sex == "Male"| dt2$Sex == "NonConclusive", ]
    Females2 <- dt1 [, c("Replicate", "AMELX", "AMELY", "Females", "Experiment")]
    Males2 <- Males2 [ , c("Replicate", "AMELX", "AMELY", "Sex", "Experiment")]
    colnames(Females2) <- c("Sample" , "AMELX", "AMELY", "Sex", "Experiment")
    colnames(Males2) <- c("Sample" , "AMELX", "AMELY", "Sex", "Experiment")
    SumTable <- rbind(Females2, Males2)
    SumTable <- SumTable[order(rank(Sample))]
    return(SumTable)
  }
  

  
########################PLOTs########################################

#####TAB 1 - STD curves#############################
  
  #Function subset the tables per experiment - (dt1 = STD2, dt2 = Reg)
  STD_Plot <- function(dt1, dt2){
    plots <- list()  # Create a list to store the individual plots
    for (i in 1:length(unique(unlist(dt1[,Experiment])))) {
      data <- dt1[dt1$Experiment == unique(unlist(dt1[,Experiment]))[i], ]
      Stat <- dt2[dt2$Experiment == unique(unlist(dt1[,Experiment]))[i], ]
      p <- plot_standard_curve(data, Stat)
      plots[[i]] <- p  # Store the plot in the list
    }
    print(do.call(grid.arrange, plots))
    
  }
  
  # Function to obtain the plot of the std dillution series per target - (dt1 = STD2, dt2 = Reg)
  plot_standard_curve <- function(dt1, dt2) {
    p <- ggplot(dt1, aes(
      x = `Normalized Area`,
      y = `Analyte Concentration`)) +
      geom_point(aes(col = `Sequence`), size = 3) +
      geom_smooth(aes(col = `Sequence`),
                  linewidth = 1,
                  method = `lm`) +
      geom_label(data = dt2,
                 position = "identity",
                 label.padding = unit(0.5, "lines"),
                 inherit.aes = FALSE,
                 aes(
                   x = -Inf,
                   y = Inf,
                   hjust = -0.2,
                   vjust = +1.1,
                   label = paste(
                     " y =",
                     formatC(Slope, digits = 2, format = "e"),
                     "x +",
                     round(Intercept, digits = 2),
                     "\n    R2=",
                     round(R2, digits = 3),
                     "\n  LOD = ",
                     formatC(LOD, digits = 2, format = "e"),
                     "\n  LOQ = ",
                     formatC(LOQ, digits = 2, format = "e")
                   )
                 ), size = 4.5) +
      facet_wrap( ~ `Sequence`, scales = "free") +
      theme_light() +
      ggtitle(paste("Standard curves")) +
      xlab("Peak area") +
      ylab("Peptide amount (pmol)") +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 35),
        strip.text = element_text(size = 18),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)
      )
    return(p)
  }
  
  #####TAB 2 - SEX IDENTIFICATION#############################
  
  #Function for the experimental based Male model - (SUMMARY = result from Sex_IDENTIFICATION, TableF = filtered data, a = outlier %)
  Plots1 <- function(SUMMARY, TableF, a){
    p1 <- ggplot(SUMMARY, aes(x= AMELX, y=AMELY))+
      geom_point(size=3, aes(col=Sex))+
      scale_color_manual(values = c("Female" = "#ffb627", "Male" = "#60d394", "NonConclusive" = "#ee6055", "Unknown" = "grey"))+
      theme_light()+
      geom_text_repel(aes(label= Sample),hjust = 0.5, vjust= -1, size= 4)+
      ggtitle("Summary")+
      theme(
        legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 30),
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 11),
        legend.text = element_text(size = 13))+
      labs(color=NULL)
    STD <- STD_Table(TableF)
    Reg <- STD_Regression(STD)
    ID_Male <- Male_ID(TableF, Reg)
    Males <- ID_Male[Sex == "Male", ]
    regressions_data2 <- Regression2(Males, a)
    FilteredMale <- as.data.table(regressions_data2$FilteredData)
    ModelInfo <-regressions_data2$ModelInfo
    p2 <- ggplot(FilteredMale, aes(x = log(AMELX, 2), y = log(AMELY, 2))) +
      geom_smooth(linewidth = 1.5, method = `lm`, col= "#EA636D", fill="#E6E7E8") +
      geom_point(size = 3, col= "#58595B") +
      geom_label(aes(
        x = -Inf,
        y = Inf,
        hjust = -0.2,
        vjust= +1.5,
        label = paste(
          " y =",
          ModelInfo$Slope,
          "x +",
          ModelInfo$Intercept,
          "\n      R2=",
          ModelInfo$R2
        ), size=4.5
      )) +
      theme_light() +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 30),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12))+
      ggtitle("Experimental Male Model")
    print(grid.arrange(p1, p2, nrow = 1, widths = c(0.6, 0.4)))
  }
  
  #Function for the predefined Male model - (fn1 = Peptide intensities data, fn2 = Experimental design, x = Choice of the modeling strategy)
  #Function for the predefined Male model - (SUMMARY = result from Sex_IDENTIFICATION)
  Plots2 <- function(SUMMARY){
    p1 <- ggplot(SUMMARY, aes(x= AMELX, y=AMELY))+
      geom_point(size=3, aes(col=Sex))+
      scale_color_manual(values = c("Female" = "#ffb627", "Male" = "#60d394", "NonConclusive" = "#ee6055", "Unknown" = "grey"))+
      theme_light()+
      geom_text_repel(aes(label= Sample),hjust = 0.5, vjust= -1, size= 4)+
      ggtitle("Summary")+
      theme(
        legend.position = "top",
        plot.title = element_text(hjust = 0.5, size = 30),
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 11),
        legend.text = element_text(size = 13))+
      labs(color=NULL)
    AMELX <- c( 24.49953985,	24.52199094,	24.70704058,	24.70972457,	25.25839351,	25.370826,	25.39675946,	25.64061101,	25.84215734,	25.88609061,
                25.95318154,	26.06919504,	26.10821049,	26.15796149,	26.23560496,	26.3146853,	26.34499341,	26.35033197,	26.36444448,	26.45871712,	
                26.56338587,	26.73711899,	26.77064635,	26.77756034,	26.79426791,	26.84062159,	26.89908305,	26.92422372,	27.0292847,	27.06167575,	
                27.09083445,	27.12415729,	27.15576323,	27.25349666,	27.274621,	27.28542113,	27.32843019,	27.42237934,	27.47421356,	27.53578819,	
                27.5829791,	27.60102059,	27.65307794,	27.96499157,	27.97749201,	27.99048029,	28.00849103,	28.13936383,	28.16542788)
    AMELY <- c( 26.14815924,	26.11695942,	26.69768664,	26.10761203,	26.5458367,	26.74468276,	27.17381855,	27.02178687,	26.9265329,	27.5351875,	
                27.47711228,	27.4777222,	27.33167096,	27.31548882,	27.50780347,	27.56764211,	27.24575753,	27.77806189,	27.62474883,	28.11807903,	
                28.00244086,	27.67873439,	28.31535925,	27.99210239,	28.16245364,	27.50159122,	27.76040322,	27.57182797,	27.8458939,	28.03370176,	
                28.30599951,	28.21038969,	27.9474296,	28.4248236,	28.16029108,	28.19554192,	28.24131775,	28.78037927,	28.69760048,	28.59320417,	
                28.43157386,	28.57224733,	28.54130137,	28.61466846,	28.90844769,	29.03409548,	28.70143946,	28.63318693,	28.7996503)
    Males <- data.table(AMELX, AMELY)
    p2 <- ggplot(Males, aes(x = log(AMELX, 2), y = log(AMELY,2))) +
      geom_smooth(linewidth = 1.5, method = `lm`, col= "#EA636D", fill="#E6E7E8") +
      geom_point(size = 3, col= "#58595B") +
      geom_label(aes(
        x = -Inf,
        y = Inf,
        hjust = -0.2,
        vjust= +1.5,
        label = paste(
          " y = 0.740x + 8.086",
          "\n     R2=0.911"
        ), size= 4.5
      )) +
      theme_light() +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 30),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12))+
      ggtitle("Predefined Male Model")
    print(grid.arrange(p1, p2, nrow = 1, widths = c(0.6, 0.4)))
    
  }
  
  ##############################MAIN FUNCTIONS###############################################
  
  #A. Plot the standard curves per experiment per target in tab 1 - (fn1 = peptide intensity data, fn2 = experimental design)
  PlotStdCurves <- function(fn1, fn2, dotProd, isoDotProd, peakRatio){
    TableF <- Load_Filter(fn1, fn2, dotProd, isoDotProd, peakRatio)
    STD <- STD_Table(TableF)
    Reg <- STD_Regression(STD)
    STD_Plot(STD, Reg) #Plot std curves
  }
  
 
  #B. Plot the overview of the AMELX and AMELY intensities per samples and the model in tab 2 - (fn1 = peptide intensity data, fn2 = experimental design, x = Choice of the modeling strategy, a = max % of values to be removed)
  Plot_Tab2 <- function(fn1, fn2, x, a, dotProd, isoDotProd, peakRatio){
    # We need to update the calls to Plots1 and Plots2 inside this function as well.
    # To do that, we'll pass the filtering parameters to the Sex_IDENTIFICATION function they call.
    SUMMARY <- Sex_IDENTIFICATION(fn1, fn2, x, a, dotProd, isoDotProd, peakRatio)
    if (x == 1) {
      # Since Plots1 calls Load_Filter again, we must update its definition and call.
      TableF <- Load_Filter(fn1, fn2, dotProd, isoDotProd, peakRatio)
      Plots1(SUMMARY, TableF, a)
    }
    else if(x == 2){
      Plots2(SUMMARY)
    }
  }
  
 
  #C. Obtain sex identification table starting from the inputs and generate a summary table in tab 2 - (fn1 = Peptide intensities data, fn2 = Experimental design, x = Choice of the modeling strategy, a = max % of values to be removed)
  Sex_IDENTIFICATION <- function(fn1, fn2, x, a, dotProd, isoDotProd, peakRatio){ 
    TableF <- Load_Filter(fn1, fn2, dotProd, isoDotProd, peakRatio)
    STD <- STD_Table(TableF) #filter the table to only have STD
    Reg <- STD_Regression(STD) #calculate LOD and LOQ per peptide and store the values in Reg
    ID_Male <- Male_ID(TableF, Reg) #Obtain male identification
    if (x == 1) { #1="Experimental_model" 
      AMELYTh <- Male_model(ID_Male, a)
    }
    else if(x == 2){ #2="Predef_model"
      AMELYTh <- Male_Predef_Model(ID_Male)
    }
    ID_Female <- Female_ID(AMELYTh, TableF) #Obtain female ID
    Summary <- Summary_table(ID_Female, ID_Male)
    return(Summary)
  }
  
  
  #D. Plot the MS traces in tab 3 -  (fn1 = Peptide intensities data, fn2 = Experimental design, fn3 = Raw intensities data, x = replicate, y = Choice of the modeling strategy, a = max % of values to be removed) 
  Plot_MSTraces <- function(fn1, fn2, fn3, x, y, a, dotProd, isoDotProd, peakRatio){ 
    SUMMARY <- Sex_IDENTIFICATION(fn1, fn2, y, a, dotProd, isoDotProd, peakRatio)
    Table <- fread(fn3)
    TableF <- Table[!grepl("precursor", Table$Transition),]
    #Filter to only keep one sample
    Plot <- TableF[Replicate == x, ]
    #Transform to long table
    Plot_Long <- Plot %>%
      separate_rows(`Raw Times`, `Raw Intensities`, sep = ",")
    #Transform to numeric
    Plot_Long$`Raw Times` <- as.numeric(Plot_Long$`Raw Times`)
    Plot_Long$`Raw Intensities` <- as.numeric(Plot_Long$`Raw Intensities`)
    #Split the table per precursor to plot the MS2 traces per peptide
    T1 <- Plot_Long[Plot_Long$`Modified Sequence` == "SIRPPYPSY", ]
    p1 <- ggplot(T1, aes(x=`Raw Times`, y=`Raw Intensities`))+
      geom_line(aes(color=Transition), size=0.6)+
      theme_light()+
      xlab("Time (min)")+
      ylab("Intensity")+
      ggtitle("AMELX-SIRPPYPSY")+
      guides(color = guide_legend(title = "Fragment"))
    
    T2 <- Plot_Long[Plot_Long$`Modified Sequence` == "SIRPPYPSYG", ]
    p2 <- ggplot(T2, aes(x=`Raw Times`, y=`Raw Intensities`))+
      geom_line(aes(color=Transition), size=0.6)+
      theme_light()+
      xlab("Time (min)")+
      ylab("Intensity")+
      ggtitle("AMELX-SIRPPYPSYG")+
      guides(color = guide_legend(title = "Fragment"))
    
    T3 <- Plot_Long[Plot_Long$`Modified Sequence` == "SM[+16]IRPPYS", ]
    p3 <- ggplot(T3, aes(x=`Raw Times`, y=`Raw Intensities`))+
      geom_line(aes(color=Transition), size=0.6)+
      theme_light()+
      xlab("Time (min)")+
      ylab("Intensity")+
      ggtitle("AMELY-SM[+16]IRPPYS")+
      guides(color = guide_legend(title = "Fragment"))
    
    T4 <- Plot_Long[Plot_Long$`Modified Sequence` == "SMIRPPY", ]
    p4 <- ggplot(T4, aes(x=`Raw Times`, y=`Raw Intensities`))+
      geom_line(aes(color=Transition), size=0.6)+
      theme_light()+
      xlab("Time (min)")+
      ylab("Intensity")+
      ggtitle("AMELY-SMIRPPY")+
      guides(color = guide_legend(title = "Fragment"))
    
    T5 <- Plot_Long[Plot_Long$`Modified Sequence` == "SM[+16]IRPPY", ]
    p5 <- ggplot(T5, aes(x=`Raw Times`, y=`Raw Intensities`))+
      geom_line(aes(color=Transition), size=0.6)+
      theme_light()+
      xlab("Time (min)")+
      ylab("Intensity")+
      ggtitle("AMELY-SM[+16]IRPPY")+
      guides(color = guide_legend(title = "Fragment"))
    
    p6 <- SUMMARY[SUMMARY$Sample == x, ]
    p6 <- as.character(p6[1,4])
    p6 <- textGrob(p6, gp = gpar(col = "black", fontsize = 40))
    
    print(grid.arrange(p6, p1, p2, p4, p5, p3, nrow = 2))
  }

  
  #########################################################
  
  
  output$samples<-renderUI({
    inFile <- input$skyline
    if (is.null(inFile))
      return(NULL)
    quant_table<-fread(inFile$datapath)
    selectInput("samples", "Select Sample", choices=levels(as.factor(quant_table$Replicate)))
  })
  
  x <- eventReactive(input$stdcurvesButton, {
    req(input$table1, input$table2)
    PlotStdCurves(input$table1$datapath, input$table2$datapath, 
                  input$dotProduct, input$isotopeDotProduct, input$peakFoundRatio)
  })
  
  v <- eventReactive(input$modelButton, {
    req(input$table1, input$table2)
    Plot_Tab2(input$table1$datapath, input$table2$datapath, 
              input$model, input$outlier,
              input$dotProduct, input$isotopeDotProduct, input$peakFoundRatio)
  })
  
  y <- eventReactive(input$sexidButton, {
    req(input$table1, input$table2)
    Sex_IDENTIFICATION(input$table1$datapath, input$table2$datapath, 
                       input$model, input$outlier,
                       input$dotProduct, input$isotopeDotProduct, input$peakFoundRatio)
  })
  
  output$download <- downloadHandler(
    filename = function() {"sex_identification_output.csv"},
    content = function(file) {
      req(input$table1, input$table2)
      result <- Sex_IDENTIFICATION(input$table1$datapath, input$table2$datapath, 
                                   input$model, input$outlier,
                                   input$dotProduct, input$isotopeDotProduct, input$peakFoundRatio)
      write.csv(result, file, row.names = FALSE)
    }
  )
  
  z <- eventReactive(input$tracesButton, {
    req(input$table1, input$table2, input$skyline)
    Plot_MSTraces(input$table1$datapath, input$table2$datapath, input$skyline$datapath, 
                  input$samples, input$model, input$outlier,
                  input$dotProduct, input$isotopeDotProduct, input$peakFoundRatio)
  })
  
  output$plot1 <- renderPlot({x()})
  output$plot2 <- renderPlot({v()})
  output$plot3 <- renderPlot({z()})
  output$result <- DT::renderDataTable({y()})
  
  })

shinyApp(ui = ui, server = server)
