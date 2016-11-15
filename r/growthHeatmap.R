#!/usr/bin/Rscript
# growthHeatmap.R
# Growth level heatmap for PMAnalyzer
#
# Author: Daniel A Cuevas
# Created on 07 May 2015
# Updated on 15 Nov 2016


# Import necessary packages
# These may need to be installed first
suppressMessages(library("ggplot2"))
suppressMessages(require("getopt"))


#################################################################
# UTILITY FUNCTIONS
#################################################################

## makeFigure()
## Create and return a ggplot object containing heatmap of growth level data
## ARGUMENTS:
##   data:           dataframe of growth curve data
##   title:          string to use as the title of the figure
##   border.color:   heatmap border color
##   low.color:      heatmap lowest gradient color
##   high.color:     heatmap highest gradient color
##   plateFlag:      boolean value specifying whether plate information
##                   is included in the dataframe
## RETURN:
##   pl:         ggplot object containing growth level heatmap
makeFigure <- function(data, title, border.color, low.color, hi.color, plateFlag) {
    # Determine growth levels for legend limits
    gl.max <- max(data$growthlevel)
    gl.min <- min(data$growthlevel)
    if (gl.min > 0.25) {
            gl.min <- 0
    }

    # Create figure
    pl <- ggplot(data, aes(x=xlab, y=sample, fill=growthlevel)) +
        geom_tile(colour=border.color) +
        theme(axis.text.y=element_text(colour="black", size=12),
              axis.text.x=element_text(colour="black", size=10,
                                       angle=90, hjust=1, vjust=0.5),
              axis.title=element_text(face="bold", size=15),
              axis.ticks=element_blank(),
              panel.grid.major=element_blank(),
              legend.key=element_rect(fill=NA),
              plot.title=element_text(face="bold")) +
        scale_x_discrete(expand=c(0, 0)) +
        scale_y_discrete(expand=c(0, 0)) +
        scale_fill_gradient(name="", low=low.color, high=hi.color,
                            breaks=c(0.25, 0.75),
                            labels=c("No Growth", "Growth"),
                            limits=c(gl.min, gl.max)) +
        ggtitle(title) + xlab("") + ylab("")

    return(pl)
}

#################################################################
# ARGUMENT PARSING
#################################################################
spec <- matrix(c(
        "infile",   "i", 1, "character",    "Input file path (required)",
        "outfile",  "o", 1, "character",    "Output file path without file type (required)",
        "plate",    "p", 0, "logical",      "Set flag if plate information is given (default: False)",
        "type",     "f", 1, "character",    "Type of image file (png*, eps, or svg) (*default)",
        "title",    "l", 1, "character",    "Title for figure (default:'')",
        "dpi",      "d", 1, "integer",      "DPI for image (default:200) (max:600)",
        "width",    "w", 1, "integer",      "Width for entire image in cm (default:45) (max:90) (required if height is specified)",
        "height",   "t", 1, "integer",      "Height for each sample in cm (default:3) (max:10) (requred if width is specified)",
        "bocol",    "b", 1, "character",    "Border color (can supply hex or name) (default:'white')",
        "locol",    "m", 1, "character",    "Low growth level color (can supply hex or name) (default:'white')",
        "hicol",    "n", 1, "character",    "High growth level color (can supply hex or name) (default:'black')",
        "help",     "h", 0, "logical",      "This help message"
        ), ncol=5, byrow=T)

opt <- getopt(spec)

# Check if help flag was given
if (!is.null(opt$help)) {
    cat(paste(getopt(spec, usage=T), "\n"))
    q(status=1)
}

# Check for input file
if (is.null(opt$infile)) {
    cat("\nInput file path not specified. Use the '-i' option.\n\n")
    cat(paste(getopt(spec, usage=T), "\n"))
    q(status=1)
}

# Check for output file
if (is.null(opt$outfile)) {
    cat("\nOutput file path not specified. Use the '-o' option.\n\n")
    cat(paste(getopt(spec, usage=T), "\n"))
    q(status=1)
}

# Check for width and height
if ((is.null(opt$width) && !is.null(opt$height))
    || (is.null(opt$height) && !is.null(opt$width))) {
    cat("\nYou must specify both height and width.\n\n")
    cat(paste(getopt(spec, usage=T), "\n"))
    q(status=1)
}

# Check image format and other metrics
if (!is.null(opt$type) && !(opt$type %in% c("png","svg","eps")) ) {
    cat("\nInvalid image type. You must specify either 'png', 'svg', or 'eps'.\n\n")
    cat(paste(getopt(spec, usage=T), "\n"))
    q(status=1)
} else if (is.null(opt$type)) {
    opt$type <- "png"
}
if (is.null(opt$dpi)) {
    opt$dpi <- 200
} else if (opt$dpi > 600) {
    opt$dpi <- 600
}
if (is.null(opt$width)) {
    opt$width <- 45  # In centimeters
} else if (opt$width > 90) {
    opt$width <- 90
}
if (is.null(opt$height)) {
    opt$height <- 3 # In centimeters
} else if (opt$height > 10) {
    opt$height <- 10
}

# Check plate flag
if (is.null(opt$plate)) {
    plateFlag <- F
} else {
    plateFlag <- opt$plate
}

# Check for title
if (is.null(opt$title)) {
    plot.title <- ""
} else {
    plot.title <- opt$title
}

# Check for border color
if (is.null(opt$bocol)) {
    border.color <- "white"
} else {
    border.color <- opt$bocol
}

# Check for low growth level color
if (is.null(opt$locol)) {
    low.color <- "white"
} else {
    low.color <- opt$locol
}

# Check for high growth level color
if (is.null(opt$hicol)) {
    hi.color <- "black"
} else {
    hi.color <- opt$hicol
}




#################################################################
# DATA PROCESSING
#################################################################
# Read in data as a table
# header=T : there is a header line
# sep="\t" : values are tab separated
# check.names=F : header names will be taken as is. There usually is a problem
#                 when numbers are part of the header
data <- read.table(opt$infile, header=T, sep="\t", check.names=F)

# Force sample names to be factor values
data$sample <- as.factor(data$sample)

# Create x-axis labels based on plate information
if (plateFlag) {
    data$xlab <- paste(data$compound, data$well, sep=" - ")
} else {
    data$xlab <- data$well
}

# Organize data by growth level
data$xlab <- factor(reorder(data$xlab, data$growthlevel, mean))

# Determine number of samples for height purposes
numS <- length(unique(data$sample))

# Add additional height to the figure if plate information is given
if (plateFlag) {
    opt$height <- numS * opt$height + 5
} else {
    opt$height <- numS * opt$height
}

pl <- makeFigure(data, plot.title, border.color, low.color, hi.color, plateFlag)
ggsave(paste(opt$outfile, ".", opt$type, sep=""),
        plot=pl,
        width=opt$width,
        height=opt$height,
        units="cm",
        dpi=opt$dpi)
