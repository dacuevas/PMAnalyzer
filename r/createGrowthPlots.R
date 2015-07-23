#!/usr/bin/Rscript
# createGrowthPlots.R
# Growth curve plot generation for PMAnalyzer
#
# Author: Daniel A Cuevas
# Created on 07 May 2015
# Updated on 07 May 2015


# Import necessary packages
# These may need to be installed first
suppressMessages(require("ggplot2"))
suppressMessages(require("scales"))
suppressMessages(require("reshape2"))
suppressMessages(require("getopt"))
suppressMessages(require("ggthemes"))
suppressMessages(require("plyr"))


#################################################################
# UTILITY FUNCTIONS
#################################################################

## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean,
## and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the
##               variable to be summariezed
##   groupvars: a vector containing names of columns
##              that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the
##                  confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {

    # New version of length which can handle NA's:
    # if na.rm==T, don't count them
	length2 <- function (x, na.rm=FALSE) {
		if (na.rm) {
            sum(!is.na(x))
        }
		else {
            length(x)
        }
	}

	# This does the summary. For each group's data frame, return a vector with
	# N, mean, and sd
    datac <- ddply(data,
                   groupvars,
                   .drop=.drop,
                   .fun = function(xx, col) {
                       c(N = length2(xx[[col]], na.rm=na.rm),
                         mean = mean(xx[[col]], na.rm=na.rm),
                         sd = sd(xx[[col]], na.rm=na.rm))
                   },
                   measurevar)

    # Rename the "mean" column
	datac <- rename(datac, c("mean" = measurevar))

    # Calculate standard error of the mean
	datac$se <- datac$sd / sqrt(datac$N)

	# Confidence interval multiplier for standard error
	# Calculate t-statistic for confidence interval:
	# e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
	ciMult <- qt(conf.interval/2 + .5, datac$N-1)
	datac$ci <- datac$se * ciMult

	return(datac)
}


## makeFigure()
## Create and return a ggplot object containing multiple growth curve plots
## ARGUMENTS:
##   plot.data:  dataframe of growth curve data
##   plot.err:   boolean value specifying whether to include standard
##               error bars in the plots
##   plateFlag:  boolean value specifying whether plate information
##               is included in the dataframe
##   title:      string to use as the title of the figure
##   color.by:   string specifying which column to color plots by
## RETURN:
##   pl:         ggplot object containing growth curve plots
makeFigure <- function(plot.data, plot.err, plateFlag, title, color.by) {
    if (plot.err) {
        grpvrs <- c("time", "well")
        if(plateFlag) {
            grpvrs <- c("mainsource", "compound", grpvrs)
        }
        # Check if a "plate" column is present
        if("plate" %in% colnames(plot.data)) {
            grpvrs <- c(grpvrs, "plate")
        }
        plot.data <- summarySE(plot.data, measurevar="od", groupvars=grpvrs)
    }

    # Create the plot using time in the x axis, optical density in the y axis.
    # Here are several aesthetic changes to the graphs
    # They include creating a line graph, a point graph, facetting based on well,
    # and other aesthetic changes to the graph with colors, labelling, etc.
    # Color is based on specified column
    if (color.by == "sample") {
        pl <- ggplot(plot.data, aes(x=time, y=od, colour=sample))
    } else if (color.by == "replicate") {
        pl <- ggplot(plot.data, aes(x=time, y=od, colour=rep))
    } else {
        pl <- ggplot(plot.data, aes(x=time, y=od))
    }
    pl <- pl + facet_wrap(~well, ncol=12)
    if (plot.err) {
        pl <- pl + geom_ribbon(aes(ymin=od-se, ymax=od+se, linetype=NA),
                                 fill="#1F77B4", alpha=0.5) +
            geom_line(size=0.8, alpha=0.5, colour="#1F77B4") +
            geom_point(size=1.5, colour="#1F77B4")
    }
    else {
        pl <- pl + geom_line(size=0.8, alpha=0.5) +
            geom_point(size=1.5, alpha=0.5)
    }
    pl <- pl + theme(axis.text=element_text(colour="black", size=12),
                     axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
                     axis.title=element_text(face="bold", size=15),
                     panel.grid.major=element_blank(),
                     panel.border=element_rect(colour="black", fill=NA),
                     strip.text=element_text(size=12),
                     strip.background=element_blank(),
                     legend.key=element_rect(fill=NA),
                     plot.title=element_text(face="bold")
        ) +
        guides(colour=guide_legend(override.aes=list(size=1)))+
    	ggtitle(title) + xlab("Time (hr)") + ylab("OD (600nm)") +
        scale_color_tableau(name="")
    return(pl)
}

#################################################################
# ARGUMENT PARSING
#################################################################
spec <- matrix(c(
        "infile",   "i", 1, "character",    "Input file path (required)",
        "outfile",  "o", 1, "character",    "Output file path without file type (required)",
        "plate",    "p", 0, "logical",      "Set flag if plate information is given (default: False)",
        "ploterr",  "e", 0, "logical",      "Set flag if error bars must be plotted (default: False)",
        "plotsep",  "s", 0, "logical",      "Set flag to create a plot per sample (default: False)",
        "colorby",  "c", 1, "character",    "Color each curve by (all*, sample, or replicate) (*default)",
        "type",     "f", 1, "character",    "Type of image file (png*, eps, or svg) (*default)",
        "title",    "l", 1, "character",    "Title for figure (appended to sample name if plotsep is specified) (default:'')",
        "dpi",      "d", 1, "integer",      "DPI for image (default:200) (max:600)",
        "width",    "w", 1, "integer",      "Width for image in cm (default:38) (max:50) (required if height is specified)",
        "height",   "t", 1, "integer",      "Height for image in cm (default:25.4) (max:50) (requred if width is specified)",
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
    opt$width <- 38  # In centimeters
} else if (opt$width > 50) {
    opt$width <- 50
}
if (is.null(opt$height)) {
    opt$height <- 25.40 # In centimeters
} else if (opt$height > 50) {
    opt$height <- 50
}

# Check plate flag
if (is.null(opt$plate)) {
    plateFlag <- F
} else {
    plateFlag <- opt$plate
}

# Check plot separation flag
if (is.null(opt$plotsep)) {
    plot.sep <- F
} else {
    plot.sep <- opt$plotsep
}

# Check which to color plots by
if (!is.null(opt$colorby) && !(opt$colorby %in% c("all", "sample", "replicate")) ) {
    cat("\nInvalid colorby type. You must specify either 'all', 'sample', or 'replicate'.\n\n")
    cat(paste(getopt(spec, usage=T), "\n"))
    q(status=1)
} else if(is.null(opt$colorby)) {
    opt$colorby <- "all"
}

# Check for error bars
if (is.null(opt$ploterr)) {
    plot.err <- F
} else {
    plot.err <- opt$ploterr
    # Also check for color by, it must be set to all
    if (opt$colorby != "all") {
        write("WARNING: Changing colorby to 'all' due to ploterr being set",
              stderr())
        opt$colorby <- "all"
    }
}


# Check for title
if (is.null(opt$title)) {
    plot.title <- ""
} else {
    plot.title <- opt$title
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

# Organize data by well numbers
data$well <- factor(data$well,
                    levels=c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11", "A12",
                             "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B11", "B12",
                             "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12",
                             "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "D11", "D12",
                             "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11", "E12",
                             "F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10", "F11", "F12",
                             "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9", "G10", "G11", "G12",
                             "H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8", "H9", "H10", "H11", "H12"))

# Change the time value from character string to numeric values
# This causes a problem when trying to plot
data$time <- as.numeric(as.character(data$time))
# Force replicate value to be a factor value
if("rep" %in% colnames(data)) {
    data$rep <- as.factor(data$rep)
}

# Create figure
if (plot.sep) {
    dataNames <- unique(data$sample)
    for (s in dataNames) {
        plot.data <- data[grep(s, data$sample, fixed=T),]
        if (plot.title == "") {
            title <- s
        } else {
            title <- paste(s, " - ", plot.title, sep="")
        }
        pl <- makeFigure(plot.data, plot.err, plateFlag, title, opt$colorby)
        ggsave(paste(opt$outfile, "_", s, ".", opt$type, sep=""),
               plot=pl,
               width=opt$width,
               height=opt$height,
               units="cm",
               dpi=opt$dpi)
    }
} else {
    pl <- makeFigure(data, plot.err, plateFlag, plot.title, opt$colorby)
    ggsave(paste(opt$outfile, ".", opt$type, sep=""),
           plot=pl,
           width=opt$width,
           height=opt$height,
           units="cm",
           dpi=opt$dpi)
}