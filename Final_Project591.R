## Author: Vamshi Mallepalli
## mvvamshi@bu.edu
## BU BF591
## Final Project 591


library(shiny)
library(ggplot2)
library(colourpicker) 
library(dplyr)
library(tidyverse)
library(bslib)               #need to install bootswatch library for theme
thematic::thematic_shiny(font = "auto")

options(shiny.maxRequestSize=30*1024^2)
js <- '.nav-tabs-custom .nav-tabs li.active {
    border-top-color: #d73925;
}"'


ui <- fluidPage(
  # Include the custom styling and theme
  theme = bs_theme(version = 4, bootswatch = "sketchy", bg = "#202123", fg = "#B8BCC2", primary = "#EA80FC", secondary = "#48DAC6"),
  titlePanel("Final Project BF591"),
  titlePanel(h3("RNA Sequence Analysis of Human Huntington Disease Brain Reveals an Extensive Increase in Inflammatory and Developmental Gene Expression App for Data Visualization", align = "center")),
  tabsetPanel(
    tabPanel("Samples",
             sidebarPanel(fileInput("input_1", paste0("Load sample data"),accept = c(".tsv",".csv")),
                          radioButtons("histo", "Chose a category to make a histogram of:", 
                                       c("PMI","age_death", "RIN","mRNA_seq_reads")),
                          colourInput("histcolor", "Histogram color"),
                          actionButton("go", "Go")
             ), #sidebarPanel
             mainPanel(
               tabsetPanel(
                 tabPanel("Summary",
                          p("This is a summary of the uploaded data."),
                          tableOutput("sample_1")
                 ),#tabPanel
                 tabPanel("Data Table",
                          p("Click"),
                          DT::dataTableOutput("sample_2")
                 ),#tabPanel
                 tabPanel("Histogram",
                          p("Here's a histogram of the data"),
                          plotOutput("sample_3")
                          
                 )
               )
             ),
    ), 
    
    tabPanel("Counts",
             sidebarPanel(
               fileInput("csv2", paste0("Load Normalized Counts")), 
               sliderInput("threshold1", label = "Choose a threshold on percentile variance", min = 0, max =100, value = 10),
               sliderInput("threshold2", label = "Choose threshold for number of samples that are non-zero", min = 0, max = 100, value = 25), 
               radioButtons("button4", "PC for the x-axis", c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")),
               radioButtons("button5", "PC for the y-axis", c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")),
             ),
             
             mainPanel(
               tabsetPanel(
                 tabPanel("Summary", textOutput("text")), 
                 tabPanel("Diagonistic Scatter Plots", plotOutput("scatter1"), plotOutput("scatter2")), 
                 tabPanel("Principal Component Analysis", plotOutput("PCscatter")),
                 tabPanel("Clustered Heatmap", plotOutput("heatmap1")))
             ),
    ),
    
    tabPanel("DE",
             sidebarPanel(
               fileInput("csv3", paste0("Load in Differential Expression Results")), 
               helpText("Here are the options for formatting the restults tab"),
               radioButtons("button8", "Chosse a column to sort the results table by", c("baseMean", "log2FoldChange", "lfcSE","stat","pvalue","padj","lfcSE")),
               helpText("Here are the options for formatting the plot and table tab"),
               helpText("A volcano plot can be generated with log2 fold-change on the x-axis and p-adjusted on the y-axis"),
               radioButtons("button", "Choose the column for the x-axis", c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")),
               radioButtons("but", "Choose the column for the y-axis", c("baseMean","log2FoldChange", "lfcSE", "stat", "pvalue", "padj")),
               colourInput("col","Base Point Color", "purple"),
               colourInput("color", "Highlight Point Color", "blue"), 
               sliderInput("range", label = "Select the magnitude of the p adjusted coloring:", min = -39, max = 0,value = -15),
               
             ),
             
             mainPanel(
               tabsetPanel(
                 tabPanel("Results", tableOutput("tabledata")),
                 tabPanel("Table", tableOutput("table2")), 
                 tabPanel("Plot", plotOutput("volcano")))
             ), 
    ),
    tabPanel("GSEA",
             sidebarPanel(
               fileInput("csv4", paste0("Load in FGSEA Results from the Differential Expression Data")),
               
               radioButtons("button6", "Choose which pathways you want to see", c("All", "positive", "negative")),
               radioButtons("button7", "Choose the column you want to sort the table by", c("pathway","pval", "padj", "log2err", "NES","size")),
               helpText("For the scatterplot tab here are the options:"),
               sliderInput("threshold5", label = "Select the adjusted p-value to filter the data", min = -.3, max = .3, value = 0.01)),
             
             mainPanel(
               tabsetPanel(
                 #abPanel("Barplot", plotOutput("barplot")),
                 tabPanel("Pathway Table", tableOutput("pathtable")),
                 tabPanel("Scatterplot", plotOutput("NEScatter")))
             ),
    ),
  ),
)

#load in the csv file, functions and output. 
server <- function(input, output, session) {
  
  #######------------------------SAMPLES TAB --------------------###########
  #load the data
  load_data4 <- reactive({
    read.csv(input$input_1$datapath,row.names = 1) #Reactive to file upload
  })
  #sample summary function
  sample_summary<- function(dataf){
    categories <- colnames(dataf) #find all the column names 
    type <- c(typeof(dataf$sample_id),typeof(dataf$PMI),typeof(dataf$age_death),
              typeof(dataf$RIN),typeof(dataf$mRNA_seq_reads)) #find the type of data stored in each column
    mean <-c("identifier",round(mean(dataf$PMI %>% na.omit),3),round(mean(dataf$age_death),3),
             round(mean(dataf$RIN),3),round(mean(dataf$mRNA_seq_reads),3)) #find the mean of each numerical column rounding to 3 digits
    sd <- ifelse(sapply(dataf, is.numeric) == TRUE,round(apply(dataf,2,sd,na.rm=TRUE),3),"identifier") #repeat with standard deviation
    summary <- data.frame(categories,type,mean,sd) #compile it into one dataframe
    return(summary)
  }
  #plot the histogram
  histogram_plot <- function(df,category,bins){
    df$category<- as.numeric(df[[category]]) #use the [[]] to change a string identifier into a column name 
    p <- ggplot(df, aes(category)) +
      geom_histogram(fill = "blue",bins = bins) #make a histogram of this data 
    return(p)
  }
  #Rendering the output
  output$text1 <- renderText({paste("Please give it a little time to process.")}) #Just a helper text
  output$sample_1 <- renderTable({req(input$input_1)  #to create the summary table output
    df<- load_data4() 
    p<- sample_summary(df)
    return(p)
  }, height = 300)
  
  output$sample_2 <- DT::renderDataTable({req(input$input_1) #to create a sortable dataframe
    DT::datatable(load_data4())})
  
  output$sample_3 <- renderPlot({req(input$input_1) #to make the histogram
    df <- load_data4()
    p<- hist(as.numeric(df[[input$histo]]), xlab= paste(input$histo), main = paste("Histogram of ",input$histo), col = input$histcolor,
             breaks=10) 
    return(p)},
    height = 300)
  
  showdata <- function(data,col_name){
    if(is.null(input$csv))
      return(NULL)
    else 
      return(showdata)
    
  }
  
  
  ##########--------------------Counts TAB---------------------############
  
  load_data2 <- reactive({
    if(is.null(input$csv2))
      return(NULL)
    df <- read.csv(input$csv2$datapath, header = TRUE, sep = ",")
    return(df)
  })
  
  #tab 2 Counts Matrix Functions 
  summarized <- function(data,threshold1, threshold2) {
    if(is.null(input$csv2))
      return(NULL)
    row.names(data) <- data$X
    data <- select(data, -X)
    data <- log(data) #log the data
    data_var <-  mutate(data, variance = apply(data,1,var)) #get the variance of the data
    data_ordered<-data_var[order(data_var$variance,decreasing = T),] #order so that you can get the top results for the percentile filtering
    thres <- threshold1/(100*length(data)) #get the threshold based on the length of the dataset 
    passed <- data_ordered[data_ordered$variance > thres,] #figure out which ones pass the threshold 
    passed <- select(passed, -variance) #drop the variance column so it doesnt mess things up with extra columns 
    data_new <- select(data_ordered, -variance)
    #filter 2 is include genes with at least X samples that are non-zero
    #so count when the value is greater than 0 and when its over the threshold keep it if not get rid of the row 
    thres2<- threshold2/(100*length(data))
    data_nonzero <- mutate(data_new, zeros = rowSums(data_var == 0))
    passed2 <- data_new[data_new$zeros > thres2,]
    total <- rbind(passed, passed2) %>% distinct() #combine the two filters to get the total and remove the duplicates
    num_samples <- ncol(data)
    num_genes <-nrow(data)
    #number and percent of genes passing current filter 
    num_passing <- nrow(total)
    percent_passing <- round((num_passing/num_genes)*100,3)
    num_failing <-num_genes - num_passing 
    percent_failing <- round((100 - percent_passing),3)
    text <- (paste0("The number of samples is ", num_samples , " and the number of genes is ", num_genes, ".", "The number of genes passing the current filter is ", num_passing, " and the percent passing is ", percent_passing,"%.", "The number of genes failing the current filter is ", num_failing , " and the percent failing is ", percent_failing, "%."))
    return(text)
  }
  
  scatter1 <- function(data, threshold1){
    if(is.null(input$csv2))
      return(NULL)
    row.names(data) <- data$X
    data <- select(data, -X)
    data <- log10(data) #follow the same process as above this is the scatter plot for the first filter with varinace
    data_var <-  mutate(data, variance = apply(data,1,var), median = apply(data,1, median))
    data_ordered<-data_var[order(data_var$variance,decreasing = T),]
    thres <- threshold1/(100*length(data))
    passed <- data_ordered[data_ordered$variance > thres,] 
    labelleddata <- mutate(data_ordered, passf1 = ifelse((row.names(data) %in% row.names(passed)), "TRUE", "FALSE")) #assign true and false to get the labels for the colors on what passed on not 
    graph <- ggplot(labelleddata, aes(x=median, y = variance, color = passf1)) +
      geom_point()+
      xlab("Log Median Count") +
      ylab("Log Variance") +
      ggtitle("Scatter Plot of Median Count vs Variance ") 
    #had alot of issues with this plot not sure if it is correct but I tried my hardest, got different result when i logged the data in the plot instead of the beggining and log scaled y 
    return(graph)
  }
  
  scatter2 <- function(data, threshold2, color){
    if(is.null(input$csv2))
      return(NULL)
    row.names(data) <- data$X
    data <- select(data, -X) #this is the second scatter plot for the second filter non zeros 
    labelleddata <- mutate(data, median = apply(data,1,median), zeros = rowSums(data == 0), passf2 = (ifelse(zeros < threshold2/(100*length(data)), "TRUE", "FALSE")))
    graph <- ggplot(labelleddata, aes(x=log(median), y = zeros , color = passf2)) +
      geom_point()+
      xlab("Median Count") +
      ylab("Number of Zeros") +
      ggtitle("Scatter Plot of Median Count vs Number of Zeros ")
    return(graph)
  }
  

  #PCA plots
  plot_pca <- function(data, x_name, y_name) {
    if(is.null(input$csv2))
      return(NULL)
    print(x_name)
    row.names(data) <- data$X
    data <- select(data, -X)
    transposed <- t(data) #transpose the data 
    pca <- prcomp(transposed)
    sd<- pca$sdev
    variance <- (sd)^2/sum((sd)^2) 
    pcs <- colnames(pca$rotation)
    cum<- cumsum(variance)
    pca_tibble <- tibble(variance_explained = variance, principal_components = as.factor(pcs), cumulative = cum)
    pca <- as.data.frame(pca$x)
    colx <- pca[,c(x_name)]
    coly <- pca[,c(y_name)]
    pc_num_x <-strtoi(substr(x_name,3,3))
    pc_num_y <- strtoi(substr(y_name,3,3))
    pcaplot <- ggplot(pca, aes(x=colx, y = coly)) +
      geom_point() +
      ggtitle("Scatter Plot of PCA Plot Projections") +
      xlab(paste0(x_name, " ", round(variance[pc_num_x]*100), "% variance")) +
      ylab(paste0(y_name, " ", round(variance[pc_num_y]*100), "% variance"))
    #need to add in the varinace for each 
    return(pcaplot)
  }
  
  #Heatmap plot 
  heatmaplot <- function(data,threshold1, threshold2){
    if(is.null(input$csv2))
      return(NULL)
    row.names(data) <- data$X  
    data <- select(data, -X)
    data <- log(data)
    data_var <-  mutate(data, variance = apply(data,1,var))
    data_ordered<-data_var[order(data_var$variance,decreasing = T),]
    thres <- threshold1/(100*length(data))
    passed <- data_ordered[data_ordered$variance > thres,]
    passed <- select(passed, -variance)
    data_new <- select(data_ordered, -variance)
    thres2<- threshold2/(100*length(data))
    data_nonzero <- mutate(data_new, zeros = rowSums(data_var == 0))
    passed2 <- data_nonzero[data_nonzero$zeros > thres2,]
    passed2 <- select(passed2, -zeros)
    total <- rbind(passed, passed2) %>% distinct()
    data1<- as.matrix(total)
    heat<- heatmap(data1, Rowv=NA, Colv=NA, col = heat.colors(256), scale="column", margins=c(5,10))  
    legend(x="right", legend=c("min", "med", "max"),fill=heat.colors(3))
    return(heat)
  }
  
  #Output functions
  output$text <- renderText(summarized(load_data2(), input$threshold1, input$threshold2))
  output$scatter1 <- renderPlot(scatter1(load_data2(), input$threshold1))
  output$scatter2 <- renderPlot(scatter2(load_data2(), input$threshold2))
  output$heatmap1 <- renderPlot(heatmaplot(load_data2(), input$threshold1, input$threshold2))
  output$PCscatter <- renderPlot(plot_pca(load_data2(), input$button4, input$button5))
  
 
  
  ############--------------------DE TAB--------------------################ 
  
  load_data1 <- reactive({
    if(is.null(input$csv3))
      return(NULL)
    df <- read.csv(input$csv3$datapath, header = TRUE, sep = ",")
    return(df)
  })
  
  draw_table1 <- function(dataf,col_name) {
    if(is.null(input$csv3))
      return(NULL)
    dataf %>% mutate(pvalue = formatC(.$pvalue, digits = 4, format = 'e'),
                     padj = formatC(.$padj, digits = 4, format = 'e')) 
    col_name <- dataf[,c(col_name)] #specify a column by name 
    sorted <- dataf[order(col_name),]
    return(sorted[1:100,])
  }
  #Volcano plot 
  volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
    if(is.null(input$csv3))
      return(NULL)
    plot<- ggplot(dataf,aes(x = !!sym(x_name), y = -log10(!!sym(y_name)))) +
      geom_point(aes(colour = cut(-log(padj), c(-Inf, -(slider), Inf)))) +
      xlab(x_name) +
      ylab(y_name) +
      ggtitle("Volcano Plot for DESeq") +
      scale_color_manual(values = c(color1, color2))
    return(plot)
  }
  
  draw_table <- function(dataf, slider) {
    if(is.null(input$csv3))
      return(NULL)
    filt_dataf <- filter(dataf, padj< 10^slider) %>%  
      mutate(pvalue = formatC(.$pvalue, digits = 4, format = 'e'),
             padj = formatC(.$padj, digits = 4, format = 'e'))
    return(filt_dataf)
  }
  
  #Output Functions
  output$volcano <-renderPlot(volcano_plot(load_data1(),input$button, input$but, input$range, input$col, input$color))
  output$table2 <- renderTable(draw_table(load_data1(), input$range))
  output$tabledata <- renderTable(draw_table1(load_data1(), input$button8))
  
  
  #############---------------------GSEA TAB-----------------###############
  
  
  load_data3 <- reactive({
    if(is.null(input$csv4))
      return(NULL)
    df <- read.csv(input$csv4$datapath, header = TRUE, sep = ",")
    return(df)
  })
  
  path_table <- function(data, threshold, but, col_name){
    if(is.null(input$csv4))
      return(NULL)
    filtered <- data[data$padj > threshold,]
    if(but == "positive"){
      fil <- filter(filtered, NES > 0)
      col_name <- fil[,c(col_name)] #specify a column by name 
      sorted <- fil[order(col_name),]
      return(sorted)
    }
    else if(but == "negative"){
      fil <- filter(filtered, NES < 0)
      col_name <- fil[,c(col_name)] #specify a column by name 
      sorted <- fil[order(col_name),]
      return(sorted)
    }
    else{
      col_name <- filtered[,c(col_name)] #specify a column by name 
      sorted <- filtered[order(col_name),]
      return(sorted)
    }
  }
  
  Scatter <- function(data, threshold){
    if(is.null(input$csv4))
      return(NULL)
    passed <- data[data$padj > threshold, "pathway"]
    newdat <- mutate(data, pass = ifelse(data$pathway %in% passed, "TRUE", "FALSE"))
    scat <- ggplot(newdat,aes(NES, -log10(padj), color = pass)) +
      geom_point()+
      scale_color_manual(values = c("#00AFBB", "grey")) +
      xlab("NES") +
      ylab("-log10(padj)") +
      ggtitle("Scatter Plot of NES vs -log10(padj) ")
    return(scat)
  }
  
  output$pathtable <- renderTable(path_table(load_data3(),input$threshold4, input$button6, input$button7))
  output$NEScatter <- renderPlot(Scatter(load_data3(),input$threshold5))
  
  
}
  
shinyApp(ui = ui, server = server)
