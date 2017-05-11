#!/usr/bin/Rscript
# box_and_density.R
# Box plots and density plots of various parameters
#
# Author: Daniel A Cuevas
# Created on 26 Sep 2016
# Updated on 14 Apr 2017

# Import necessary packages
# These may need to be installed first
suppressMessages(require("getopt"))
suppressMessages(require("ggplot2"))
suppressMessages(require("ggthemes"))
suppressMessages(require("reshape2"))
suppressMessages(require("grid"))
suppressMessages(require("plyr"))


#################################################################
# UTILITY FUNCTIONS
#################################################################
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

makePlots <- function(data, title=NULL, ftitle=NULL) {
    if (is.null(title)) {
        title <- "All Samples"
        ftitle <- "all_samples"
    }
    # Calculate statistics
    data.stats <- summarySE(data, measurevar="val", groupvars=c("metric"))

    # Draw box plots
    pl <- ggplot(data, aes(x=x, y=val)) +
        facet_wrap(~metric, ncol=4, scales="free") +
        geom_boxplot() +
        theme(axis.text=element_text(colour="black", size=10),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title=element_text(face="bold", size=15),
              panel.grid.minor=element_blank(),
              panel.border=element_rect(colour="black", fill=NA),
              panel.spacing=unit(3, "mm"),
              strip.text=element_text(face="bold", size=12, vjust=0),
              strip.background=element_rect(colour="white",
                                            fill=NA, size=3),
              legend.key=element_rect(fill=NA),
              plot.title=element_text(face="bold")) +
        scale_y_continuous(labels=fmt_decimals(3)) +
        xlab("") + ylab("") + ggtitle(title)

    ggsave(paste(opt$outpath, "/box_plots_", ftitle,".png", sep=""),
           plot=pl,
           width=20,
           height=20,
           units="cm",
           dpi=200)

    # Draw density plots
    pl <- ggplot(data, aes(x=val)) +
        facet_wrap(~metric, ncol=4, scales="free") +
        geom_density(adjust=0.25, fill="#c7c7c7") +
        geom_rect(data=data.stats, aes(x=NULL, y=NULL, xmin=val-sd, xmax=val+sd,
                                       ymin=-Inf, ymax=Inf), alpha=0.2, fill="#1F77B4") +
        geom_vline(aes(xintercept=val), data.stats, colour="#1F77B4", size=1) +
        theme(axis.text=element_text(colour="black", size=10),
              axis.title=element_text(face="bold", size=15),
              panel.grid.minor=element_blank(),
              panel.border=element_rect(colour="black", fill=NA),
              panel.spacing=unit(3, "mm"),
              strip.text=element_text(face="bold", size=12, vjust=0),
              strip.background=element_rect(colour="white",
                                            fill=NA, size=3),
              legend.key=element_rect(fill=NA),
              plot.title=element_text(face="bold")) +
        scale_x_continuous(labels=fmt_decimals(3)) +
        xlab("") + ylab("") + ggtitle(title)

    ggsave(paste(opt$outpath, "/density_plots_", ftitle ,".png", sep=""),
           plot=pl,
           width=35,
           height=20,
           units="cm",
           dpi=200)
}

# Taken from http://stackoverflow.com/questions/10035786/ggplot2-y-axis-label-decimal-precision
fmt_decimals <- function(decimals=0) {
    function(x) format(x, digits=decimals, nsmall=3)
}

#################################################################
# ARGUMENT PARSING
#################################################################
spec <- matrix(c(
        "infile",     "i", 1, "character",    "Input file path (required)",
        "outpath",    "o", 1, "character",    "Output file path (required)",
        "plate",      "p", 0, "logical",      "Set flag if plate information is given (default: False)",
        "help",       "h", 0, "logical",      "This help message"
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

# Check for output filepath
if (is.null(opt$outpath)) {
    cat("\nOutput filepath not specified. Use the '-o' option.\n\n")
    cat(paste(getopt(spec, usage=T), "\n"))
    q(status=1)
}

# Check for output suffix
if (is.null(opt$out_suffix)) {
    out.suffix <- "out"
} else {
    out.suffix <- opt$out_suffix
}

# Check plate flag
if (is.null(opt$plate)) {
    plateFlag <- F
} else {
    plateFlag <- opt$plate
}

#################################################################
# DATA PROCESSING
#################################################################
data <- read.delim(opt$infile)

# Force sample and replicate to be factor values
data$sample <- as.factor(data$sample)
data$rep<- as.factor(data$rep)

if (plateFlag) {
    idvars <- c("sample", "rep", "well", "mainsource", "compound")
} else {
    idvars <- c("sample", "rep", "well")
}
melt.data <- melt(subset(data, select=c(-growthclass)),
                  id.vars=idvars, variable.name=c("metric"), value.name="val")
melt.data$x <- 0

# Make plots of all data
makePlots(melt.data)

# Make plots for each sample
for (s in unique(melt.data$sample)) {
    subdata <- subset(melt.data, sample == s)
    makePlots(subdata, s, s)
}
