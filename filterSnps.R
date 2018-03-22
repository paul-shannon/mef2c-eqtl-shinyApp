library(GenomicRanges)
library(RUnit)
library(MEF2C.data)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("mef2c"))
   mef2c <- MEF2C.data()
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_shoulder.0()
   test_shoulder.0.rosmap()
   test_get.tf.bindingSites()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
intersectTables <- function(tbl.foreground, tbl.background)
{
   stopifnot(colnames(tbl.foreground)[1:3] == c("chrom", "start", "end"))
   stopifnot(colnames(tbl.background)[1:3] == c("chrom", "start", "end"))
   tbl.foreground <- tbl.foreground[order(tbl.foreground$start, decreasing=FALSE),]
   tbl.background <- tbl.background[order(tbl.background$start, decreasing=FALSE),]

   tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl.foreground), GRanges(tbl.background)))
   colnames(tbl.ov) <- c("foreground", "background")

   tbl.foreground[unique(tbl.ov$foreground),]

} # intersectTables
#------------------------------------------------------------------------------------------------------------------------
applyFilters <- function(session, filters, trena.model, snpShoulder)
{
   stopifnot(trena.model %in% c("model_cer", "model_tcx", "model_ros"))
   stopifnot(snpShoulder >= 0 & snpShoulder < 250)
   stopifnot(all(filters %in% c("footprints", "dhs", "enhancers")))

   tbl.activeSnps <-  mef2c@misc.data[["IGAP.snpChip"]] # 228
   tbl.activeSnps$start <- tbl.activeSnps$start - snpShoulder
   tbl.activeSnps$end <- tbl.activeSnps$start + snpShoulder

   tbl.enhancers <- mef2c@misc.data$enhancer.locs
   study.roi <- list(chrom="chr5", start=min(tbl.enhancers$start) - 10000, end=max(tbl.enhancers$end) + 10000)
   load("tbl.fp.geneSymbolsAndShortMotifAdded.RData")

   tbl.activeFp <- data.frame()
   tbl.activeFpInModel <- data.frame()

   #if("footprints" %in% filters){
   tbl.activeSnps <- intersectTables(tbl.activeSnps, tbl.fp)
   tbl.activeFp <- intersectTables(tbl.fp, tbl.activeSnps)
    # }

   tbl.model <- switch(trena.model,
      "model_cer" = getModels(mef2c)[["mef2c.cory.wgs.cer.tfClass"]],
      "model_tcx" = getModels(mef2c)[["mef2c.cory.wgs.tcx.tfClass"]],
      "model_ros" = getModels(mef2c)[["mef2c.cory.wgs.ros.tfClass"]])

   tfoi <- tbl.model$gene[1:10]

   if(nrow(tbl.activeFp) > 0)
      tbl.activeFpInModel <- subset(tbl.activeFp, geneSymbol.tfclass %in% tfoi | geneSymbol.motifdb %in% tfoi)

   rownames(tbl.activeFpInModel) <- NULL

   tbl.activeSnpsInModel <- data.frame()
   if(nrow(tbl.activeSnps) > 0 & nrow(tbl.activeFpInModel) > 0){
      tbl.activeSnpsInModel <- intersectTables(tbl.activeSnps, tbl.activeFpInModel)
      }

   tbl.activeFpInModel$name <- with(tbl.activeFpInModel,
                          sprintf("%s.%s.%s", shortMotif, geneSymbol.tfclass, geneSymbol.motifdb))

   tbl.cleanedUp <- tbl.activeFpInModel[, c("chrom", "start", "end", "name", "score", "strand")]

   if(!is.null(session)){
      temp.filename <- sprintf("igvdata/tmp%d.bed", as.integer(Sys.time()))
      write.table(tbl.cleanedUp, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
      model.name <- sub("model_", "", trena.model)
      trackName <- sprintf("snp+%d/fp/%s", snpShoulder, model.name)
      session$sendCustomMessage(type="displayBedTrack",
                                message=(list(filename=temp.filename,
                                              trackName=trackName,
                                              color="blue",
                                              trackHeight=40)))
      } # is session

   invisible(list(fp=tbl.activeFpInModel, snp=tbl.activeSnpsInModel, fpClean=tbl.cleanedUp))

} # applyFilters
#------------------------------------------------------------------------------------------------------------------------
test_shoulder.0 <- function()
{
   printf("--- test_shoulder.0")
   x <- applyFilters(session=NULL, filters="footprints", trena.model="model_cer", snpShoulder=0)
   checkEquals(names(x), c("fp", "snp", "fpClean"))
   xyz <- 99

} # test_shoulder.0
#------------------------------------------------------------------------------------------------------------------------
test_shoulder.0.rosmap <- function()
{
   printf("--- test_shoulder.0.rosmap")
   x <- applyFilters(session=NULL, filters="footprints", trena.model="model_ros", snpShoulder=0)
   checkEquals(names(x), c("fp", "snp", "fpClean"))
   xyz <- 99

} # test_shoulder.0
#------------------------------------------------------------------------------------------------------------------------
get.tf.bindingSites <- function(tf, roi=NULL)
{
   if(!exists("tbl.fpMappedToAllTfs"))
      load("tbl.fpMappedToAllTFs.RData")    # put this in MEF2C.data

   tbl.tf.bindingSites <- subset(tbl.fpMappedToAllTFs, geneSymbol.motifdb == tf | geneSymbol.tfclass==tf)
   tbl.tfbs <- tbl.tf.bindingSites
   if(!is.null(roi)){
      tbl.roi <- as.data.frame(roi, stringsAsFactors=FALSE)
      tbl.tfbs <- intersectTables(tbl.tf.bindingSites, tbl.roi)
      }
   tbl.tfbs <- tbl.tfbs[, c("chrom", "start", "end", "shortMotif")]
   tbl.tfbs <- tbl.tfbs[order(tbl.tfbs$start, decreasing=FALSE),]
   rownames(tbl.tfbs) <- NULL
   tbl.tfbs

} # get.tf.bindingSites
#------------------------------------------------------------------------------------------------------------------------
test_get.tf.bindingSites <- function()
{
   printf("--- test_get.tf.bindingSite.for.model")

   tbl.tfbs <- get.tf.bindingSites("EGR3", roi=list(chrom="chr5", start=88883402, end=88883979))
   checkEquals(dim(tbl.tfbs), c(21, 4))
   checkEquals(length(unique(tbl.tfbs$shortMotif)), 10)

      # now the whole region covered by our tbl.fpMappedToAllTFs

   tbl.tfbs <- get.tf.bindingSites("EGR3")
   checkEquals(dim(tbl.tfbs), c(1216, 4))
   checkEquals(length(unique(tbl.tfbs$shortMotif)), 18)

} # test_get.tf.bindingSite.for.model
#------------------------------------------------------------------------------------------------------------------------
# EGR3 is the 4th-ranked MEF2C tf in model #2, mayo temporal cortex
# one of the snps (IGAP?  Mayo?) at 5:88,883,759-88,883,758,  within an EGR3 binding site about  1k downstream
# of the somewhat non-sensical biomart-reported mef2c promoter
# see screenshot in my email to cory (17 mar 2018) titled "two screenshots from the shiny app"
# reproduce this here, ensuring that footprint, model, and snps are all in sync
#
# > x
# $fp
#   chrom    start      end                          motifName strand    score length distance.from.tss                                                            id      db geneSymbol.orig pubmedID shortMotifName shortMotif geneSymbol.motifdb pubmedID.1 geneSymbol.tfclass pubmedID.2                name
# 1  chr5 88692233 88692242  Hsapiens-jaspar2016-KLF5-MA0599.1      + 10.37080     10           -192233  MEF2C.fp.downstream.192233.Hsapiens-jaspar2016-KLF5-MA0599.1 hint_16            KLF5 24194598       MA0599.1   MA0599.1               KLF5   24194598               EGR3   23180794  MA0599.1.EGR3.KLF5
# 2  chr5 88692233 88692243   Hsapiens-jaspar2016-SP1-MA0079.3      + 10.61800     11           -192233   MEF2C.fp.downstream.192233.Hsapiens-jaspar2016-SP1-MA0079.3 hint_16             SP1 24194598       MA0079.3   MA0079.3                SP1   24194598               EGR3   23180794   MA0079.3.EGR3.SP1
# 3  chr5 88883756 88883773 Hsapiens-jaspar2016-KLF13-MA0657.1      -  7.35955     18              -710 MEF2C.fp.downstream.000710.Hsapiens-jaspar2016-KLF13-MA0657.1 well_16           KLF13 24194598       MA0657.1   MA0657.1              KLF13   24194598               EGR3   23180794 MA0657.1.EGR3.KLF13
# 4  chr5 88883757 88883771   Hsapiens-jaspar2016-SP2-MA0516.1      - 17.14610     15              -709   MEF2C.fp.downstream.000709.Hsapiens-jaspar2016-SP2-MA0516.1 well_16             SP2 24194598       MA0516.1   MA0516.1                SP2   24194598               EGR3   23180794   MA0516.1.EGR3.SP2
#
# $snp
#         chrom    start      end         id    score mut ref    beta stderr     pval  hg19pos
# 2362907  chr5 88692240 88692240  rs6894196 2.241921   T   C  0.1711 0.0619 0.005729 87988057
# 2363210  chr5 88883758 88883758 rs80043958 2.616903   G   A -0.0622 0.0205 0.002416 88179575
#
# $fpClean
#   chrom    start      end                name    score strand
# 1  chr5 88692233 88692242  MA0599.1.EGR3.KLF5 10.37080      +
# 2  chr5 88692233 88692243   MA0079.3.EGR3.SP1 10.61800      +
# 3  chr5 88883756 88883773 MA0657.1.EGR3.KLF13  7.35955      -
# 4  chr5 88883757 88883771   MA0516.1.EGR3.SP2 17.14610      -
test_snpsIntersectingModels <- function()
{
   printf("--- test_snpsIntersectingModels")
   x <- applyFilters(session=NULL, filters="footprints", trena.model="model_tcx", snpShoulder=0)

   tbl.motifs <- mef2c@misc.data[["allDNA-jaspar2018-human-mouse-motifs"]]
   tfs <- getModels(mef2c)[["mef2c.cory.wgs.tcx.tfClass"]]$gene[1:5]
   for(tf in tfs){
     motifs.this.tf <- geneToMotif(MotifDb, tf, "tfclass")$motif
     coi <- c("chrom", "start", "end", "name", "score")
     tbl.motifs.tf <- subset(tbl.motifs, name %in% motifs.this.tf)[, coi]
     tbl.motifs.tf <- tbl.motifs.tf[order(tbl.motifs.tf$start, decreasing=FALSE),]
     temp.filename <- sprintf("igvdata/tmp%d.bed", as.integer(Sys.time()))
     write.table(tbl.motifs.tf, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, file=temp.filename)
     }


} # test_snpsIntersectingModesl
#------------------------------------------------------------------------------------------------------------------------
