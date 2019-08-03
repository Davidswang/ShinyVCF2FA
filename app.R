#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(seqinr)
library(stringr)

#Define UI functions for file input tools
VCF_UI <- function(id) {
    ns = NS(id)
    
    list(
        textOutput(ns("output_area")), 
        fileInput("file1", "Choose VCF File",
                  multiple = FALSE,
                  accept = ".vcf")
    )
}

fasta_UI <- function(id) {
    ns = NS(id)
    
    list(
        textOutput(ns("output_area")), 
        fileInput("file2", "Choose Reference FASTA File",
                  multiple = FALSE,
                  accept = c(".fa",".fasta"))
    )
}

# Define UI for VCF to FASTA App
# Sidebar panel for user interaction (upload/download)
# Main panel to preview data/results

ui<-fluidPage(
    # App title ----
    titlePanel("VCF to FASTA Processing"),
    
    # Sidebar layout with input and output definitions ----
    sidebarLayout(
        
        # Sidebar panel for inputs ----
        sidebarPanel(
            
            # Input: Select a file ----
            VCF_UI("test1"),
            actionButton("pvButton", "Preview vcf file"),
            # Horizontal line, provides spacing
            tags$hr(),
            # Input: Select a file ----
            fasta_UI("test2"),
            # Horizontal line, provides spacing
            tags$b("Download modified FASTA file\n"),
            tags$hr(),
            downloadButton('downloadData', 'Download')
        ),
        
        # Main panel for displaying outputs ----
        mainPanel(
            
            # Output: Preview Original FASTA 
            h4("FASTA sequence sliced from first variant to last variant"),
            div(style="width:500px;",verbatimTextOutput("value2")),
            
            # Output: Preview Modified FASTA
            h4("Alternative FASTA sequence, variants capitalized"),
            div(style="width:500px;",verbatimTextOutput("value3")),
            
            # Output: Preview VCF file
            h4("Preview Data from VCF File"),
            tableOutput("contents")
        )
    )
) 


# Define server logic required to draw a histogram
server <- function(input, output) {

    output$contents <- renderTable({
        
        # input$file1 will be NULL initially. After the user selects
        # and uploads a file, head of that data file by default,
        # or all rows if selected, will be shown.
        
        req(input$file1)
        req(input$pvButton)
        
        # when reading semicolon separated files,
        # having a comma separator causes `read.csv` to error
        tryCatch(
            {
                df <- read.table(input$file1$datapath,
                                 header = FALSE,
                                 sep = '\t')
                names(df)<- c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','VALUES')
            },
            error = function(e) {
                # return a safeError if a parsing error occurs
                stop(safeError(e))
            }
        )
        #return(head(df))
        return(df)
    })
    
    output$table <- renderTable({
        datasetInput()
    })
    
    output$value <- renderText({
        
        # input$file1 will be NULL initially. After the user selects
        # and uploads a file, head of that data file by default,
        # or all rows if selected, will be shown.
        
        req(input$file1)
        
        # when reading semicolon separated files,
        # having a comma separator causes `read.csv` to error
        tryCatch(
            {
                df <- read.table(input$file1$datapath,
                                 header = FALSE,
                                 sep = '\t')
                names(df)<- c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','VALUES')
            },
            error = function(e) {
                # return a safeError if a parsing error occurs
                stop(safeError(e))
            }
        )
        
        startBase<-df[1,2]
        last <- nrow(df)
        lastBase<-df[last,2]
        
        myString <- paste('First variant called at', startBase,'and last variant called at',lastBase)
        
        return(myString)
    })
    
    output$value2 <- renderText({
        req(input$file1)
        req(input$file2)
        
        tryCatch(
            {
                df <- read.table(input$file1$datapath,
                                 header = FALSE,
                                 sep = '\t')
                names(df)<- c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','VALUES')
            },
            error = function(e) {
                # return a safeError if a parsing error occurs
                stop(safeError(e))
            }
        )
        
        startBase<-df[1,2]
        last <- nrow(df)
        lastBase<-df[last,2]
        
        tryCatch(
            {
                df2 <- read.fasta(input$file2$datapath, as.string=TRUE)
            },
            error = function(e) {
                # return a safeError if a parsing error occurs
                stop(safeError(e))
            }
        )
        
        test<-df2[[1]]
        myString <-gsub(" ","",test, fixed=FALSE)
        myString2 <-substr(myString, startBase,lastBase)
        #myString <-"TESTING"
        return(myString2)
    })
    
    output$value3 <- renderText({
        req(input$file1)
        req(input$file2)
        
        tryCatch(
            {
                df2 <- read.fasta(input$file2$datapath, as.string=TRUE)
            },
            error = function(e) {
                # return a safeError if a parsing error occurs
                stop(safeError(e))
            }
        )
        
        tryCatch(
            {
                df <- read.table(input$file1$datapath,
                                 header = FALSE,
                                 sep = '\t')
                names(df)<- c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','VALUES')
            },
            error = function(e) {
                # return a safeError if a parsing error occurs
                stop(safeError(e))
            }
        )
        
        startBase<-df[1,2]
        last <- nrow(df)
        lastBase<-df[last,2]
        
        #myString <- paste('File uploaded for',names(df2)[[1]])
        # myString <-substr(df2[[1]],startBase,lastBase)
        # 
        # #altBase <- df[1,5]
        
        
        test<-df2[[1]]
        myString <-gsub(" ","",test, fixed=FALSE)
        myString2 <-substr(myString, startBase,lastBase)
        #substr(myString2,1,1) <- as.character(df[1,5])
        #myString <-"TESTING"
        #return(myString2)
        
        for (i in 1:last){
            pos<-df[i,2]-startBase+1
            substr(myString2, pos, pos) <- as.character(df[i,5])
        }
        
        return(myString2)
        
        
    })
    
    thedata <- renderText({
        req(input$file1)
        req(input$file2)
        
        tryCatch(
            {
                df2 <- read.fasta(input$file2$datapath, as.string=TRUE)
            },
            error = function(e) {
                # return a safeError if a parsing error occurs
                stop(safeError(e))
            }
        )
        
        tryCatch(
            {
                df <- read.table(input$file1$datapath,
                                 header = FALSE,
                                 sep = '\t')
                names(df)<- c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','VALUES')
            },
            error = function(e) {
                # return a safeError if a parsing error occurs
                stop(safeError(e))
            }
        )
        
        startBase<-df[1,2]
        last <- nrow(df)
        lastBase<-df[last,2]
        
        #myString <- paste('File uploaded for',names(df2)[[1]])
        # myString <-substr(df2[[1]],startBase,lastBase)
        # 
        # #altBase <- df[1,5]
        
        
        test<-df2[[1]]
        myString <-gsub(" ","",test, fixed=FALSE)
        myString2 <-substr(myString, startBase,lastBase)
        
        for (i in 1:last){
            pos<-df[i,2]-startBase+1
            substr(myString2, pos, pos) <- as.character(df[i,5])
        }
        
        #myString3 <- as.dataframe(myString2)
        #names(myString3)<-">Modified_FASTA"
        return(myString2)
        
        
    })
    
    # downloadHandler() takes two arguments, both functions.
    # The content function is passed a filename as an argument, and
    #   it should write out data to that filename.
    output$downloadData <- downloadHandler(
        
        # This function returns a string which tells the client
        # browser what name to use when saving the file.
        filename = function() {
            paste("Modified", "fa", sep = ".")
        },
        
        # This function should write data to a file given to it by
        # the argument 'file'.
        content = function(file) {
            # Write to a file specified by the 'file' argument
            write.table(thedata(), file, sep = "\t",col.names=">Modified_FASTA",row.names = FALSE,quote=FALSE)
        }
    )

}

# Run the application 
shinyApp(ui = ui, server = server)
