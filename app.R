library(shiny)
library(DT)
library(igvShiny)
library(htmlwidgets)
library(GenomicRanges)
#----------------------------------------------------------------------------------------------------
library(MEF2C.data)

if(!exists("mef2c"))
   mef2c <- MEF2C.data()

source("filterSnps.R")

addResourcePath("igvdata", "igvdata") # so the shiny webserver can see, and return track files for igv
#----------------------------------------------------------------------------------------------------
currentTracks <- list()
currentFilters <- list()
currentShoulder <- 0
currentModel <- "model_cer"
filter.result <- list(fp=data.frame(), snp=data.frame(), fpClean=data.frame())
#----------------------------------------------------------------------------------------------------
ui <- fluidPage(

  includeScript("message-handler.js"),
  tags$head(tags$link(rel = "stylesheet", type = "text/css",
                     href = "http://maxcdn.bootstrapcdn.com/font-awesome/4.2.0/css/font-awesome.min.css")
     ),

  titlePanel("Evaluate IGAP eQTLs"),

  sidebarLayout(
     sidebarPanel(
        width=3,
        selectInput("targetGene", "Target gene:",
                    c("MEF2C" = "MEF2C",
                      "TREM2" = "TREM2")),

        radioButtons("roi", "Genomic Region of Interest:",
                     c(#"Set Manually" = "manualSelection",
                       "TSS +/- 2kb" = "tss_2kb",
                       "TSS +/- 5kb" = "tss_5kb",
                       "Study region" = "studyRegion")),

        checkboxGroupInput("tracks", "Genome Tracks - toggle display",
                           c("Enhancers" = "enhancers",
                             "TSS" = "tss",
                             "DHS" = "dhs",
                             "Footprints" = "fp",
                             "Top 5 TF binding sites model 1" = "tfbs.model.1",
                             "Top 5 TF binding sites model 2" = "tfbs.model.2",
                             "Top 5 TF binding sites model 3" = "tfbs.model.3"
                             )),

        #radioButtons("motifTfMapping", "Map Motifs to TFs Using",
        #             c("MotifDb (parsimonious)" = "motifdb",
        #               "TFClass (expansive)" = "tfclass",
        #               "Both" = "both")),

        #checkboxGroupInput("filters", "Filters",
        #                   c("Enhancers" = "enhancers",
        #                     "Footprints" = "footprints")),

        radioButtons("model.trn", "Choose Model",
                    c("Model 1 (Mayo cerebellum)" = "model_cer",
                      "Model 2 (Mayo temporal cortex)" = "model_tcx",
                      "Model 3 (ROSMAP)" = "model_ros")),


        sliderInput("snpShoulder",
                    "SNP shoulder",
                    value = 0,
                    min = 0,
                    max = 100),
        actionButton("showSNPsButton", "Show all SNPs"),
        actionButton("findSnpsInModelButton", "Find intersecting SNPs")
        ),

     mainPanel(
       tabsetPanel(type="tabs",
                   id="trenaTabs",
                   tabPanel(title="IGV",     value="igvTab",       igvShinyOutput('igvShiny')),
                   tabPanel(title="Summary", value="summaryTab",   verbatimTextOutput("summary")),
                   tabPanel(title="TRN",    value="snpTableTab",  DTOutput("snpSummaryTable")),
                   tabPanel(title="debug",   value="debugTab",     verbatimTextOutput("debug"))
                   )
       ) # mainPanel
    ) # sidebarLayout
  ) # fluidPage

#--------------------------------------------------------------------------------
server <- function(input, output, session) {

    # define a reactive conductor. it returns a function, is
    # redefined with a change to dist or n
    # has reactive children  plot, summary, table, below, which mention it
    # thereby establishing the reactive chain.

   observeEvent(input$showSNPsButton, {
      temp.filename <- calculateVisibleSNPs(isolate(input$targetGene),
                                     isolate(input$roi),
                                     isolate(input$filters),
                                     isolate(input$snpShoulder))
     session$sendCustomMessage(type="displaySnps", message=(list(filename=temp.filename)))
     })

   observeEvent(input$findSnpsInModelButton, {
      trena.model <- isolate(input$model.trn)
      filters <- isolate(input$filters)
      snpShoulder <- isolate(input$snpShoulder)
      filter.result <- applyFilters(session, filters, trena.model, snpShoulder)
      })

   observeEvent(input$snpShoulder, {
      currentShoulder <- input$snpShoulder
      })

   observeEvent(input$tracks, {
     # updateTabsetPanel(session, "trenaTabs", selected="debugTab")
     trackNames <- isolate(input$tracks)
     printf("currentTracks: %s", paste(currentTracks, collapse=","))
     newTracks <- setdiff(trackNames, currentTracks)
     printf("newTracks: %s", paste(newTracks, collapse=","))
     if(length(newTracks) > 0){
        output$debug <- renderPrint({newTracks})
        displayTrack(session, newTracks)
        }
     currentTracks <<- sort(unique(unlist(c(currentTracks, newTracks))))
     })

   observeEvent(input$filters, {
      currentFilters <- input$filters
      print(paste(currentFilters, collapse=","))
      })

   observeEvent(input$model.trn, {
      currentModel <- input$model.trn
      print(paste(currentFilters, collapse=","))
      })

   observeEvent(input$roi, {
     currentValue = input$roi
     tss <- 88884466;   # get from MEF2C.data
     roi.string <- switch(currentValue,
                          tss_2kb = sprintf("chr5:%d-%d", tss-2000, tss+2000),
                          tss_5kb = sprintf("chr5:%d-%d", tss-5000, tss+5000),
                          studyRegion = "chr5:88391000-89322000")
     session$sendCustomMessage(type="roi", message=(list(roi=roi.string)))
     })


  output$summary <- renderPrint({
    print(input$targetGene)
    print(input$roi)
    print(input$filters)
    print(input$snpShoulder)
    print(input$model.trn)
    })

   output$snpSummaryTable <- renderDT({
     tbl <- filter.result$fp #mtcars
     model.name <- input$model.trn
     tbl.model <- switch(model.name,
        "model_cer" = getModels(mef2c)[["mef2c.cory.wgs.cer.tfClass"]],
        "model_tcx" = getModels(mef2c)[["mef2c.cory.wgs.tcx.tfClass"]],
        "model_ros" = getModels(mef2c)[["mef2c.cory.wgs.ros.tfClass"]])

     tbl.model <- cbind(cbind(tbl.model[, 1], as.data.frame(lapply(tbl.model[, 2:10], function(col) signif(col, 3)))), tbl.model[, 11:12])
     colnames(tbl.model)[1] <- "TF"
     return(tbl.model)
     })

   output$igvShiny <- renderIgvShiny({
     igvShiny("hello shinyApp")
     })

  } # server

#--------------------------------------------------------------------------------
calculateVisibleSNPs <- function(targetGene, roi, filters, snpShoulder)
{
   #return(sprintf("calculateVisibleSNPs(%s), %s, %s, %d", date(), targetGene, roi, snpShoulder))
   tbl.mayo.snps <- mef2c@misc.data$MAYO.eqtl.snps[, c("chrom", "start", "end", "score")]
   tbl.mayo.snps <- tbl.mayo.snps[order(tbl.mayo.snps$start, decreasing=FALSE),]

   tbl.igap.snps <- mef2c@misc.data[["IGAP.snpChip"]][, c("chrom", "start", "end", "score")]
   tbl.igap.snps <- tbl.igap.snps[order(tbl.igap.snps$start, decreasing=FALSE),]

   temp.filename <- sprintf("igvdata/tmp%d.bedGraph", as.integer(Sys.time()))
   write.table(tbl.igap.snps, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
   return(temp.filename)

} # calculateVisibleSNPs
#--------------------------------------------------------------------------------
displayTrack <- function(session, trackNames)
{
   if("enhancers" %in% trackNames){
      tbl.enhancers <- mef2c@misc.data$enhancer.locs
      temp.filename <- sprintf("igvdata/tmp%d.bed", as.integer(Sys.time()))
      write.table(tbl.enhancers, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
      session$sendCustomMessage(type="displayBedTrack",
                                message=(list(filename=temp.filename,
                                              trackName="enhancers",
                                              color="black",
                                              trackHeight=40)))
      } # enhancers

   if("tss" %in% trackNames) {
      tss <- mef2c@misc.data$TSS
      tbl.tss <- data.frame(chrom="chr5", start=tss, end=tss, stringsAsFactors=FALSE)
      temp.filename <- sprintf("igvdata/tmp%d.bed", as.integer(Sys.time()))
      write.table(tbl.tss, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
      session$sendCustomMessage(type="displayBedTrack",
                                message=(list(filename=temp.filename,
                                              trackName="TSS",
                                              color="blue",
                                              trackHeight=40)))
     } # tss

    if("dhs" %in% trackNames){
      tbl.dhs <- mef2c@misc.data$tbl.dhs
      temp.filename <- sprintf("igvdata/tmp%d.bed", as.integer(Sys.time()))
      write.table(tbl.dhs, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
      session$sendCustomMessage(type="displayBedTrack",
                                message=(list(filename=temp.filename,
                                              trackName="DHS",
                                              color="darkGreen",
                                              trackHeight=40)))
      } # dhs

    if("fp" %in% trackNames){
      roi <- list(chrom="chr5", start=88391000, end=89322000)
      tbl.fp <- getFootprints(mef2c, roi)[, c("chrom", "start", "end", "shortMotifName", "score")]
      temp.filename <- sprintf("igvdata/tmp%d.bed", as.integer(Sys.time()))
      write.table(tbl.fp, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
      session$sendCustomMessage(type="displayBedTrack",
                                message=(list(filename=temp.filename,
                                              trackName="Footprints",
                                              color="gray",
                                              trackHeight=40)))
      } # fp

    if("tfbs.model.1" %in% trackNames){
      tbl.motifs <- mef2c@misc.data[["allDNA-jaspar2018-human-mouse-motifs"]]
      tfs <- getModels(mef2c)[["mef2c.cory.wgs.cer.tfClass"]]$gene[1:5]
      for(tf in tfs){
        tbl.tfbs <- get.tf.bindingSites(tf)
        temp.filename <- sprintf("igvdata/tmp%d.bed", as.integer(Sys.time()))
        write.table(tbl.tfbs, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
        session$sendCustomMessage(type="displayBedTrack",
                                message=(list(filename=temp.filename,
                                              trackName=tf,
                                              color="darkRed",
                                              trackHeight=40)))
        } # for tf
      } # tfbs.model.1

    if("tfbs.model.2" %in% trackNames){
      tbl.motifs <- mef2c@misc.data[["allDNA-jaspar2018-human-mouse-motifs"]]
      tfs <- getModels(mef2c)[["mef2c.cory.wgs.tcx.tfClass"]]$gene[1:5]
      for(tf in tfs){
        tbl.tfbs <- get.tf.bindingSites(tf)
        temp.filename <- sprintf("igvdata/tmp%d.bed", as.integer(Sys.time()))
        write.table(tbl.tfbs, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
        session$sendCustomMessage(type="displayBedTrack",
                                message=(list(filename=temp.filename,
                                              trackName=tf,
                                              color="darkRed",
                                              trackHeight=40)))
        } # for tf
      } # tfbs.model.2

    if("tfbs.model.3" %in% trackNames){
      tbl.motifs <- mef2c@misc.data[["allDNA-jaspar2018-human-mouse-motifs"]]
      tfs <- getModels(mef2c)[["mef2c.cory.wgs.ros.tfClass"]]$gene[1:5]
      for(tf in tfs){
        tbl.tfbs <- get.tf.bindingSites(tf)
        temp.filename <- sprintf("igvdata/tmp%d.bed", as.integer(Sys.time()))
        write.table(tbl.tfbs, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
        session$sendCustomMessage(type="displayBedTrack",
                                message=(list(filename=temp.filename,
                                              trackName=tf,
                                              color="darkRed",
                                              trackHeight=40)))
        } # for tf
      } # tfbs.model.3

} # displayTrack
#------------------------------------------------------------------------------------------------------------------------
# shinyApp(ui, server)
