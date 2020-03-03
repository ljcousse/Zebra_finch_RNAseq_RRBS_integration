#################
# RCircos plots # 
#################


# Load in Libraries
library(RCircos)
library(xlsx)
library(biomaRt)

# Load in chromosome info
df <- read.table("genome.txt",stringsAsFactors = F,header=F)

# Load in gene annotation for Zebra finch genome
ensembl=useMart("ensembl")
ensembl = useDataset("tguttata_gene_ensembl",mart=ensembl)
genes <- getBM(attributes = c('chromosome_name','start_position','end_position','hgnc_symbol','ensembl_gene_id'),mart = ensembl)
print(dim(genes))
genes <- genes[!duplicated(genes$ensembl_gene_id),]
print(dim(genes))

# Test data in RCircos package
# data(UCSC.HG19.Human.CytoBandIdeogram);
# head(UCSC.HG19.Human.CytoBandIdeogram);


# Create cyto info object
cyto.info <- cbind(df[,1],rep(0,nrow(df)),df[,2],rep("",nrow(df)),rep("black",nrow(df)))
colnames(cyto.info) <- c("Chromosome","ChromStart","ChromEnd","Band","Stain")
cyto.info <- as.data.frame(cyto.info,stringsAsFactors=F)
cyto.info$ChromStart <- as.integer(cyto.info$ChromStart)
cyto.info$ChromEnd <- as.integer(cyto.info$ChromEnd)

# Set core Components
chr.exclude <- NULL
tracks.inside <- 5
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)

# Get plot parameters 
rcircos.params <- RCircos.Get.Plot.Parameters();
rcircos.cyto <- RCircos.Get.Plot.Ideogram();
rcircos.position <- RCircos.Get.Plot.Positions();

# Set some adjustmetns
rcircos.params$heatmap.width <- 10
RCircos.Reset.Plot.Parameters(rcircos.params);
RCircos.List.Plot.Parameters()

# Load in DE genes
age_genes <- read.table("GeneLists/ListGenes_age.txt",header=F,stringsAsFactors = F)$V1
age_genes_df <- genes[genes$ensembl_gene_id%in%age_genes,c(1:3,5)]
age_genes_

sex_genes <- unique(c(read.table("GeneLists/ListGenes_1dph.txt",header=F,stringsAsFactors = F)$V1,
                      read.table("GeneLists/ListGenes_20dph.txt",header=F,stringsAsFactors = F)$V1,
                      read.table("GeneLists/ListGenes_65dph.txt",header=F,stringsAsFactors = F)$V1,
                      read.table("GeneLists/ListGenes_adult.txt",header=F,stringsAsFactors = F)$V1))
sex_genes_df <- genes[genes$ensembl_gene_id%in%sex_genes,c(1:3,5)]
sex_genes_df <- sex_genes_df[!grepl("RRCB",sex_genes_df$chromosome_name),]

#Get df for genes in the right form
age_genes_df$chromosome_name <- toupper(age_genes_df$chromosome_name)
age_genes_df$chromosome_name <- paste("chr",age_genes_df$chromosome_name,sep="")
age_genes_df <- age_genes_df[!grepl("RANDOM",age_genes_df$chromosome_name),]
age_genes_df <- age_genes_df[!grepl("UN",age_genes_df$chromosome_name),]
age_genes_df <- age_genes_df[!grepl("RRCB",age_genes_df$chromosome_name),]
# age_genes_df <- age_genes_df[!grepl("LG5",age_genes_df$chromosome_name),]
# age_genes_df <- age_genes_df[!grepl("LG2",age_genes_df$chromosome_name),]
age_genes_df$linedata <- rep(1,nrow(age_genes_df))
colnames(age_genes_df) <- c("Chromosome","chromStart","chromEnd","EnsGeneName","LineData")


sex_genes_df$chromosome_name <- toupper(sex_genes_df$chromosome_name)
sex_genes_df$chromosome_name <- paste("chr",sex_genes_df$chromosome_name,sep="")
sex_genes_df <- sex_genes_df[!grepl("RANDOM",sex_genes_df$chromosome_name),]
sex_genes_df <- sex_genes_df[!grepl("UN",sex_genes_df$chromosome_name),]
sex_genes_df$linedata <- rep(2,nrow(sex_genes_df))
colnames(sex_genes_df) <- c("Chromosome","chromStart","chromEnd","EnsGeneName","LineData")



##############################################
# Rewritten fucntions of the RCircos package #
##############################################

# Rewrite to get different color for chromosomal bars
#######################################################

# Set color you want as background for your chromosome plot
color_bg <- "coral3"



RCircos.Draw.Chromosome.Ideogram_L <- function (ideo.pos = NULL, ideo.width = NULL) 
{
  RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
  RCircos.Pos <- RCircos.Get.Plot.Positions()
  RCircos.Par <- RCircos.Get.Plot.Parameters()
  if (is.null(ideo.pos)) 
    ideo.pos <- RCircos.Par$chr.ideo.pos
  if (is.null(ideo.width)) 
    ideo.width <- RCircos.Par$chrom.width
  outerPos <- ideo.pos + ideo.width
  innerPos <- ideo.pos
  chromosomes <- unique(RCircos.Cyto$Chromosome)
  RCircos.Track.Outline(outerPos, innerPos, num.layers = 1, 
                        chromosomes, track.colors = rep(color_bg, length(chromosomes)))
  whiteBands <- which(RCircos.Cyto$BandColor == "white")
  darkBands <- RCircos.Cyto[-whiteBands, ]
  for (aBand in seq_len(nrow(darkBands))) {
    aColor <- darkBands$BandColor[aBand]
    aStart <- darkBands$StartPoint[aBand]
    aEnd <- darkBands$EndPoint[aBand]
    posX <- c(RCircos.Pos[aStart:aEnd, 1] * outerPos, RCircos.Pos[aEnd:aStart, 
                                                                  1] * innerPos)
    posY <- c(RCircos.Pos[aStart:aEnd, 2] * outerPos, RCircos.Pos[aEnd:aStart, 
                                                                  2] * innerPos)
    polygon(posX, posY, col = aColor, border = NA)
  }
}

RCircos.Chromosome.Ideogram.Plot_L <- function (tick.interval = 0) 
{
  RCircos.Draw.Chromosome.Ideogram_L()
  RCircos.Highligh.Chromosome.Ideogram()
  if (tick.interval > 0) {
    RCircos.Ideogram.Tick.Plot(tick.interval)
  }
  RCircos.Label.Chromosome.Names()
}



# Rewrite to go outside standard heatmap colors
################################################

RCircos.Get.Heatmap.Color.Scale_L <- function (heatmap.color = NULL) 
{
  if (is.null(heatmap.color)) 
    stop("Missing argument in RCircos.Get.Heatmap.Color.Scale().\n")
  allOnes <- seq(1, 1, length = 256)
  allZeros <- seq(0, 0, length = 256)
  one2zeor <- seq(1, 0, length = 256)
  zero2one <- seq(0, 1, length = 256)
  if (heatmap.color == "BlueWhiteRed") {
    RedRamp <- rgb(allOnes, one2zeor, one2zeor)
    BlueRamp <- rgb(zero2one, zero2one, allOnes)
    ColorRamp <- cbind(BlueRamp, RedRamp)
  }
  else if (heatmap.color == "GreenWhiteRed") {
    RedRamp <- rgb(allOnes, one2zeor, one2zeor)
    GreenRamp <- rgb(zero2one, allOnes, zero2one)
    ColorRamp <- cbind(GreenRamp, RedRamp)
  }
  else if (heatmap.color == "GreenYellowRed") {
    RedRamp <- rgb(allOnes, one2zeor, allZeros)
    GreenRamp <- rgb(zero2one, allOnes, allZeros)
    ColorRamp <- cbind(GreenRamp, RedRamp)
  }
  else if (heatmap.color == "GreenBlackRed") {
    RedRamp <- rgb(zero2one, allZeros, allZeros)
    GreenRamp <- rgb(allZeros, one2zeor, allZeros)
    ColorRamp <- cbind(GreenRamp, RedRamp)
  }
  else if (heatmap.color == "YellowToRed") {
    ColorRamp <- rgb(allOnes, one2zeor, allZeros)
  } 
  else if (heatmap.color == "goldenrod1") {
    ColorRamp <- rep("goldenrod1",256)
  }
  else if (heatmap.color == "goldenrod2") {
    ColorRamp <- rep("goldenrod2",256)
  }
  else if (heatmap.color == "goldenrod3") {
    ColorRamp <- rep("goldenrod3",256)
  }
  else {
    ColorRamp <- rgb(one2zeor, one2zeor, one2zeor)
  }
  return(ColorRamp)
}



RCircos.Heatmap.Plot_L <- function (heatmap.data = NULL, data.col = NULL, track.num = NULL, 
                                    side = c("in", "out"), min.value = NULL, max.value = NULL, 
                                    inside.pos = NULL, outside.pos = NULL, genomic.columns = 3, 
                                    is.sorted = TRUE) 
{
  if (is.null(heatmap.data)) 
    stop("Genomic data missing in RCircos.Heatmap.Plot().\n")
  if (is.null(genomic.columns)) 
    stop("Missing number of columns for genomic position.\n")
  if (is.null(data.col) || data.col <= genomic.columns) 
    stop("Data column must be ", genomic.columns + 1, " or bigger.\n")
  boundary <- RCircos.Get.Plot.Boundary(track.num, side, inside.pos, 
                                        outside.pos, FALSE)
  outerPos <- boundary[1]
  innerPos <- boundary[2]
  RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
  RCircos.Pos <- RCircos.Get.Plot.Positions()
  RCircos.Par <- RCircos.Get.Plot.Parameters()
  colorMap <- RCircos.Get.Heatmap.Color.Scale_L(RCircos.Par$heatmap.color)
  if (is.null(min.value) || is.null(max.value)) {
    columns <- (genomic.columns + 2):ncol(heatmap.data)
    min.value <- min(as.matrix(heatmap.data[, columns]))
    max.value <- max(as.matrix(heatmap.data[, columns]))
  }
  colorLevel <- seq(min.value, max.value, length = length(colorMap))
  heatmap.data <- RCircos.Get.Single.Point.Positions(heatmap.data, 
                                                     genomic.columns)
  plotLocations <- RCircos.Get.Start.End.Locations(heatmap.data, 
                                                   RCircos.Par$heatmap.width)
  chromosomes <- unique(as.character(RCircos.Cyto$Chromosome))
  outlineColors <- rep("white", length(chromosomes))
  RCircos.Track.Outline(outerPos, innerPos, num.layers = 1, 
                        chrom.list = chromosomes, track.colors = outlineColors)
  heatmapValues <- as.numeric(heatmap.data[, data.col])
  for (aPoint in 1:length(heatmapValues)) {
    theLevel <- which(colorLevel >= heatmapValues[aPoint])
    cellColor <- colorMap[min(theLevel)]
    theStart <- plotLocations[aPoint, 1]
    theEnd <- plotLocations[aPoint, 2]
    polygonX <- c(RCircos.Pos[theStart:theEnd, 1] * outerPos, 
                  RCircos.Pos[theEnd:theStart, 1] * innerPos)
    polygonY <- c(RCircos.Pos[theStart:theEnd, 2] * outerPos, 
                  RCircos.Pos[theEnd:theStart, 2] * innerPos)
    polygon(polygonX, polygonY, col = cellColor, border = NA)
  }
}






RCircos.Reset.Plot.Parameters_L <- function (new.params = NULL) 
{
  if (is.null(new.params)) 
    stop("Missing function argument.\n")
  old.params <- RCircos.Get.Plot.Parameters()
  if (new.params$radius.len != old.params$radius.len || new.params$plot.radius != 
      old.params$plot.radius || new.params$chr.ideo.pos != 
      old.params$chr.ideo.pos || new.params$tracks.inside != 
      old.params$tracks.inside || new.params$tracks.outside != 
      old.params$tracks.outside) {
    stop("Please use RCircos.Set.Core.Components() instead.\n")
  }
  if (new.params$chr.ideo.pos != old.params$chr.ideo.pos || 
      new.params$highlight.pos != old.params$highlight.pos || 
      new.params$chr.name.pos != old.params$chr.name.pos) {
    stop("Please use customized ideogram plot methods instead.\n")
  }
  RCircos.Validate.Plot.Parameters_L(new.params)
  if (new.params$chrom.width != old.params$chrom.width) {
    differ <- new.params$chrom.width - old.params$chrom.width
    new.params$highlight.pos <- old.params$highlight.pos + 
      differ
    new.name.pos <- old.params$chr.name.pos + differ
    if (new.params$chr.name.pos < new.name.pos) 
      new.params$chr.name.pos <- new.name.pos
  }
  if (new.params$track.in.start >= new.params$chr.ideo.pos) 
    new.params$track.in.start <- new.params$chr.ideo.pos - 
    0.05
  new.name.end <- new.params$chr.name.pos + 0.3
  if (new.params$track.out.start < new.name.end) 
    new.params$track.out.start <- new.name.end
  if (new.params$track.padding != old.params$track.padding || 
      new.params$track.height != old.params$track.height) {
    message(paste0("Track height and/or track padding have been ", 
                   "reset\n. Actual total data track plotted may differ.\n"))
  }
  if (old.params$base.per.unit != new.params$base.per.unit && 
      old.params$chrom.paddings == new.params$chrom.paddings) {
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
    band.len <- RCircos.Cyto$ChromEnd - RCircos.Cyto$ChromStart
    genome.len <- sum(as.numeric(band.len))
    padding.const <- RCircos.Get.Padding.Constant()
    total.units <- genome.len/new.params$base.per.unit
    new.padding <- round(padding.const * total.units, digits = 0)
    if (new.padding != new.params$base.per.unit) {
      message(paste("\nNote: chrom.padding", new.params$chrom.paddings, 
                    " was reset to", new.padding, "\n"))
      new.params$chrom.paddings <- new.padding
    }
  }
  RCircosEnvironment <- NULL
  RCircosEnvironment <- get("RCircos.Env", envir = globalenv())
  RCircosEnvironment[["RCircos.PlotPar"]] <- NULL
  RCircosEnvironment[["RCircos.PlotPar"]] <- new.params
  if (old.params$base.per.unit != new.params$base.per.unit || 
      old.params$chrom.paddings != new.params$chrom.paddings) {
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
    RCircos.Cyto <- RCircos.Cyto[, 1:5]
    RCircosEnvironment[["RCircos.Cytoband"]] <- NULL
    RCircos.Set.Cytoband.Data(RCircos.Cyto)
    RCircosEnvironment[["RCircos.Base.Position"]] <- NULL
    RCircos.Set.Base.Plot.Positions()
  }
}



RCircos.Validate.Plot.Parameters_L <- function (parameters = NULL) 
{
  if (is.null(parameters)) 
    stop("Missing function argument.\n")
  for (aParam in seq_along(parameters)) {
    if (class(parameters[[aParam]]) == "numeric" && parameters[[aParam]] < 
        0) {
      stop("Plot parameters cannot have negative values.")
    }
  }
  # heatmapColors <- RCircos.Get.Supported.HeatmapColors()
  # if ((parameters$heatmap.color %in% heatmapColors) == FALSE) {
  #   stop("\nUnsupported heatmap color defined.\nType \n            RCircos.Get.Supported.HeatmapColors() for more info.")
  # }
  colorNames <- colors()
  if (!parameters$text.color %in% colorNames || !parameters$hist.color %in% 
      colorNames || !parameters$line.color %in% colorNames || 
      !parameters$scatter.color %in% colorNames || !parameters$tile.color %in% 
      colorNames) {
    stop("Unsupported R plot color defined.")
  }
  background <- parameters$track.background
  if (!is.null(background) && !is.na(background)) {
    if (!background %in% colorNames) 
      stop("Unsupported color for track background.")
  }
  grid.color <- parameters$grid.line.color
  if (!is.null(grid.color) && !is.na(grid.color)) {
    if (!grid.color %in% colorNames) 
      stop("Unsupported color for grid line.")
  }
}





RCircos.Chromosome.Ideogram.Plot_jpg <- function (tick.interval = 0) 
{
  RCircos.Draw.Chromosome.Ideogram_L()
  RCircos.Highligh.Chromosome.Ideogram()
  if (tick.interval > 0) {
    RCircos.Ideogram.Tick.Plot(tick.interval)
  }
  RCircos.Label.Chromosome.Names_jpg()
}



RCircos.Label.Chromosome.Names_jpg <- function (chr.name.pos = NULL) 
{
  RCircos.Par <- RCircos.Get.Plot.Parameters()
  if (is.null(chr.name.pos)) {
    chr.name.pos <- RCircos.Par$chr.name.pos
  }
  else {
    if (chr.name.pos < RCircos.Par$track.in.start) 
      stop("Chromosome name positions overlap with inside track.\n")
    if (chr.name.pos > RCircos.Par$track.out.start) 
      message("May not plot data tracks outside of chromosomes.\n")
  }
  RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
  RCircos.Pos <- RCircos.Get.Plot.Positions()
  chroms <- unique(RCircos.Cyto$Chromosome)
  rightSide <- nrow(RCircos.Pos)/2
  for (aChr in seq_len(length(chroms))) {
    chrRows <- which(RCircos.Cyto$Chromosome == chroms[aChr])
    chrStart <- min(RCircos.Cyto$StartPoint[chrRows])
    chrEnd <- max(RCircos.Cyto$EndPoint[chrRows])
    mid <- round((chrEnd - chrStart + 1)/2, digits = 0) + 
      chrStart
    chrName <- sub(pattern = "chr", replacement = "", chroms[aChr])
    text(RCircos.Pos[mid, 1] * RCircos.Par$chr.name.pos, 
         RCircos.Pos[mid, 2] * RCircos.Par$chr.name.pos, label = chrName, 
         pos = ifelse(mid <= rightSide, 4, 2), srt = 0,cex=3)
  }
}


####################
# Production plots #
####################

# Storage as pdf
####################

# Open pdf
out.file <- "RCircos_Plastoscine_DEgenes.pdf";
pdf(file=out.file, height=8, width=8, compress=TRUE);

# Remove boxes 
par(lty="blank")

# Set plot area
RCircos.Set.Plot.Area();

# Plot choromosomes
RCircos.Chromosome.Ideogram.Plot_L()

# Plot DE genes 
rcircos.params$heatmap.color <- "goldenrod1"
RCircos.Reset.Plot.Parameters_L(rcircos.params);
RCircos.Heatmap.Plot_L(age_genes_df,track.num=1,"in",data.col=5)
rcircos.params$heatmap.color <- "goldenrod2"
RCircos.Reset.Plot.Parameters_L(rcircos.params);
RCircos.Heatmap.Plot_L(sex_genes_df,track.num=2,"in",data.col=5)
# rcircos.params$heatmap.color <- "goldenrod3"
# RCircos.Reset.Plot.Parameters_L(rcircos.params);
# RCircos.Heatmap.Plot_L(int_genes_df,track.num=3,"in",data.col=5)


dev.off()

# Storage as Jpeg
#####################

# Open jpg
out.file <- "RCircos_Plastoscine_DEgenes.jpg";
jpeg(file=out.file, height=750, width=750);

# Adjust text size for jpg
rcircos.params$text.size <- 10
rcircos.params$char.width <- 1000
RCircos.Reset.Plot.Parameters_L(rcircos.params);

# Remove boxes 
par(lty="blank")

# Set plot area
RCircos.Set.Plot.Area();

# Plot choromosomes
RCircos.Chromosome.Ideogram.Plot_L()

# Plot DE genes 
rcircos.params$heatmap.color <- "goldenrod1"
RCircos.Reset.Plot.Parameters_L(rcircos.params);
RCircos.Heatmap.Plot_L(age_genes_df,track.num=1,"in",data.col=5)
rcircos.params$heatmap.color <- "goldenrod2"
RCircos.Reset.Plot.Parameters_L(rcircos.params);
RCircos.Heatmap.Plot_L(sex_genes_df,track.num=2,"in",data.col=5)
# rcircos.params$heatmap.color <- "goldenrod3"
# RCircos.Reset.Plot.Parameters_L(rcircos.params);
# RCircos.Heatmap.Plot_L(int_genes_df,track.num=3,"in",data.col=5)


dev.off()



# #########################
# # Possible usefull code #
# #########################
# 
# # RCircos.Area.Highlight(highlight.area=c("chr1",12700000,34600000), track.num=c(1:10),
# # side=c("in"), hightlight.color=rgb(0.5, 0.5, 0, 0.5),
# # inside.pos=NULL, outside.pos=NULL)
# 
# RCircos.Gene.Connector.Plot(genes_df,track.num=1, side="in");
# RCircos.Gene.Name.Plot(genes_df,track.num = 2,side="in",name.col = 4)
# RCircos.Line.Plot(genes_df,data.col=5,1,"in")
# RCircos.Point.Plot(genes_df,track.num=4,"in",data.col=5,min.value = 0,max.value = 1000,outside.pos = 1.29,inside.pos = 1.17
#                    ,point.type=20,point.scale=1)
# 
# library(circlize)
# 
# # Create data
# data = data.frame(
#   factor = sample(letters[1:8], 1000, replace = TRUE),
#   x = rnorm(1000), 
#   y = runif(1000)
# )
# 
# # Step1: Initialise the chart giving factor and x-axis.
# circos.initialize( factors=data$factor, x=data$x )
# 
# # Step 2: Build the regions. 
# circos.trackPlotRegion(factors = data$factor, y = data$y, panel.fun = function(x, y) {
#   circos.axis()
# })
# 
# # library
# library(circlize)
# 
# # Create data
# set.seed(123)
# data = data.frame(
#   factor = sample(letters[1:8], 1000, replace = TRUE),
#   x = rnorm(1000), 
#   y = runif(1000)
# )
# 
# # Initialize the plot.
# par(mar = c(1, 1, 1, 1) ) 
# circos.initialize(factors = data$factor, x = data$x )
# 
# # Build the regions of track #1
# circos.trackPlotRegion(factors = data$factor, y=data$y , bg.col = rgb(0.1,0.1,seq(0,1,0.1),0.4) , bg.border = NA)
# 
# # Add a link between a point and another
# circos.link("a", 0, "b", 0, h = 0.4)
# 
# # Add a link between a point and a zone
# circos.link("e", 0, "g", c(-1,1), col = "green", lwd = 2, lty = 2, border="black" )
# 
# # Add a link between a zone and another
# circos.link("c", c(-0.5, 0.5), "d", c(-0.5,0.5), col = "red", border = "blue", h = 0.2)
# 
# #Create data
# name=c(3,10,10,3,6,7,8,3,6,1,2,2,6,10,2,3,3,10,4,5,9,10)
# feature=paste("feature ", c(1,1,2,2,2,2,2,3,3,3,3,3,3,3,4,4,4,4,5,5,5,5) , sep="")
# dat <- data.frame(name,feature)
# dat <- with(dat, table(name, feature))
# 
# # Charge the circlize library
# library(circlize)
# 
# # Make the circular plot
# chordDiagram(as.data.frame(dat), transparency = 0.5)
# 
# data(RCircos.Gene.Label.Data);
# name.col <- 4;
# side <- "in";
# track.num <- 1;
# RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data[RCircos.Gene.Label.Data$Chromosome!="chrX",],track.num, side);
# track.num <- 2;
# RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data,name.col,track.num, side);
