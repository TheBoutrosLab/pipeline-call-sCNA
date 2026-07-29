### Required packages
required.packages <- c(
    'argparse',
    'BoutrosLab.plotting.general',
    'naturalsort'
    );
### Check and install missing packages
missing.packages <- required.packages[!(required.packages %in% installed.packages())];
if (length(missing.packages)) {
    install.packages(missing.packages);
    }

### Load required packages and hide printing of loaded package names
invisible(
    lapply(
        X = required.packages,
        FUN = library,
        character.only = TRUE
        )
    );

pdf(NULL);

parser <- ArgumentParser();
parser$add_argument(
    '-s',
    '--sample-id',
    type = 'character',
    help = 'Sample ID of the sample type, either normal or tumor sample ID',
    required = TRUE
    );
parser$add_argument(
    '-t',
    '--tumor-id',
    type = 'character',
    help = 'Tumor ID to compare column header in BAF and LogR data',
    required = TRUE
    );
parser$add_argument(
    '-b',
    '--baf-file',
    type = 'character',
    help = 'Path to the B Allele Frequency (BAF) file of sample type normal or tumor',
    required = TRUE
    );
parser$add_argument(
    '-l',
    '--logr-file',
    type = 'character',
    help = 'Path to the Log R Ratio file of sample type normal or tumor',
    required = TRUE
    );
parser$add_argument(
    '-p',
    '--position-scale',
    type = 'character',
    choices = c('index', 'genome-position'),
    default = 'genome-position',
    help = 'Chromosome axis style, whether data points are ordered serially "index" or spatially ordered based on "genome-position"'
    );
parser$add_argument(
    '-c',
    '--per-chrom-mode',
    action = 'store_true',
    help = 'Create LogR/BAF plots per chromosome'
    );
parser$add_argument(
    '-o',
    '--output-dir',
    type = 'character',
    help = 'Output directory path to publish logR/BAF plots',
    required = TRUE
    );
parser$add_argument(
    '-r',
    '--plot-resolution',
    type = 'integer',
    help = 'Resolution of the output logR/BAF plot',
    required = TRUE
    );
parser$add_argument(
    '-d',
    '--reference-dict-file',
    type = 'character',
    help = 'Path to the reference FASTA index file',
    required = TRUE
    );

### Constants
abline.pos.baf <- 0.5;
abline.pos.logr <- NULL;
ylab.baf <- 'BAF';
ylab.logr <- expression(bold('log')[bold('2')] * bold('R'));
ylab.position.baf <- seq(0, 1, by = 0.2);
ylab.position.logr <- seq(-2, 2, by = 1);
ylim.baf <- c(-0.05, 1.05);
ylim.logr <- c(-2.2, 2.2);

### Functions ###
basepairs.to.kilobasepairs <- function(x.basepairs) {
    x.basepairs / 1000;
    }
modify.positions <- function(df, max.end.pos) {
    # get the first pos to add (max last pos in chr 1)
    pos.to.add <- max.end.pos[1 == max.end.pos$Chromosome.serial.num,'Length_kbp'];
    grid.pos <- pos.to.add;
    df[1 == df$Chromosome.serial.num,'Position.mod'] <- df[1 == df$Chromosome.serial.num,'Position'];

    for (chr in unique(df$Chromosome.serial.num)[-1]) {
        df[chr == df$Chromosome.serial.num,'Position.mod'] <- pos.to.add + df[chr == df$Chromosome.serial.num,'Position'];
        pos.to.add <- pos.to.add + max.end.pos[max.end.pos$Chromosome.serial.num == chr,'Length_kbp'];
        grid.pos <- c(grid.pos, pos.to.add);
        }
    return(list(df.mod.pos = df, grid.pos = grid.pos));
    }

generate.xlabel.position <- function(xgrid.position) {
    total.points <- length(xgrid.position) - 1;
    xlabel.position <- lapply(
        1:total.points, function(i) {
                (xgrid.position[i] + xgrid.position[i + 1]) / 2;
                }
            );
    xlabel.position <- as.numeric(xlabel.position[1:total.points]);
    return(xlabel.position);
    }

scatterplot.wrapper <- function(data, x, y, xaxis.label, xlab.position, xgrid.position, xlab, ylab, ylab.position, ylim, abline.pos)  {
    create.scatterplot(
        data = data,
        formula = get(y) ~ x,
        cex = 0.1,
        col = 'red',
        xat = xlab.position,
        xaxis.lab = xaxis.label,
        xaxis.cex = 1,
        xlab.cex = 2.5,
        xlab.label = xlab,
        xlimits = c(0, max(xgrid.position)),
        yat = ylab.position,
        yaxis.cex = 1.4,
        ylab.cex = 2.5,
        ylimits = ylim,
        ylab.label = ylab,
        ygrid.at = NULL,
        xgrid.at = xgrid.position,
        add.grid = TRUE,
        grid.colour = 'black',
        line.lwd = 2,
        abline.h = abline.pos,
        abline.col = 'white',
        resolution = plot.resolution
        );
    }

generate.multipanelplot <- function(plot.file.path, plot.title) {
    create.multipanelplot(
        filename = plot.file.path,
        plot.objects = list(
            logr.plot,
            baf.plot
            ),
        main = plot.title,
        layout.width = 1,
        layout.height = 2,
        plot.objects.heights = c(5, 5),
        y.spacing = -2,
        width = 20,
        height = 10,
        resolution = plot.resolution
        );
    }

### Parse arguments
args <- parser$parse_args();

sample.id <- args$sample_id;
tumor.id <- args$tumor_id;
variant.position.scale <- args$position_scale;
per.chrom.mode <- args$per_chrom_mode;
plot.resolution <- args$plot_resolution;
baf.file <- read.delim(
    file = args$baf_file,
    header = TRUE,
    as.is = TRUE,
    sep = '\t',
    check.names = FALSE
    );
logr.file <- read.delim(
    file = args$logr_file,
    header = TRUE,
    as.is = TRUE,
    sep = '\t',
    check.names = FALSE
    );

### Update column names
baf.name <- paste0(sample.id, '.BAF');
names(baf.file)[tumor.id == names(baf.file)] <- baf.name;
logr.name <- paste0(sample.id, '.logR');
names(logr.file)[tumor.id == names(logr.file)] <- logr.name;

### Merge BAF and logR data
baf.logr <- merge(
    x = baf.file,
    y = logr.file,
    by = c('Chromosome', 'Position')
    );

### Create df with unique chromosomes and index number
### chrX and chrY will be assigned '23' and '24' if present
chr.unique <- naturalsort(unique(baf.logr$Chromosome));
chr.serial.num <- data.frame(
    Chromosome.serial.num = 1:length(chr.unique),
    Chromosome = chr.unique
    );

### Merge chrosomsome index number with BAF/logR data
baf.logr <- merge(
    x = chr.serial.num,
    y = baf.logr,
    by = 'Chromosome'
    );

baf.logr <- baf.logr[order(baf.logr$Chromosome.serial.num, baf.logr$Position),];

### Get Chromosome Lengths
if (isTRUE(per.chrom.mode) || ('genome-position' == variant.position.scale)) {
    ref.dict <- read.table(
        file = args$reference_dict_file,
        as.is = TRUE,
        sep = '\t',
        fill = TRUE
        );
    chr.col <- grep('SN:', ref.dict);
    chr.length.col <- grep('LN:', ref.dict);
    chr.length.list <- c();
    for (chr.n in 1:length(chr.unique)) {
        chr.pattern <- paste0('SN:', chr.unique[chr.n]);
        chr.length.list[chr.n] <- ref.dict[ref.dict[,chr.col] == chr.pattern, chr.length.col];
        }
    chr.length.list <- as.numeric(sub('LN:', '', chr.length.list));
    chr.length.list <- basepairs.to.kilobasepairs(x.basepairs = chr.length.list);
    chr.length <- data.frame(
        Chromosome = chr.unique,
        Length_kbp = chr.length.list
        );
    chr.length <- merge(
        x = chr.serial.num,
        y = chr.length,
        by = 'Chromosome'
        );
    chr.length <- chr.length[order(chr.length$Chromosome.serial.num, chr.length$Length_kbp),];
    }

### Per Chromosome Dataset - Plot BAF and logR based on variant position in each chromosome
if (isTRUE(per.chrom.mode)) {
    print('Generating LogR/BAF plot per chromosome!');
    baf.logr$Position <- basepairs.to.kilobasepairs(x.basepairs = baf.logr$Position);
    for (chr in chr.unique) {
        main.title <- paste(sample.id, '-', chr);
        plot.file.path <- file.path(
            args$output_dir,
            paste0(sample.id, '_LogR-BAF-', chr, '.png')
            );

        baf.logr.per.chr <- baf.logr[chr == baf.logr$Chromosome,];

        xgrid.position <- c(0, chr.length[chr == chr.length$Chromosome, 'Length_kbp']);

        baf.plot <- scatterplot.wrapper(
            data = baf.logr.per.chr,
            x = baf.logr.per.chr$Position,
            y = baf.name,
            xaxis.label = NA,
            xlab.position = TRUE,
            xgrid.position = xgrid.position,
            ylab.position = ylab.position.baf,
            xlab = 'Position (kbp)',
            ylab = ylab.baf,
            ylim = ylim.baf,
            abline.pos = abline.pos.baf
            );

        logr.plot <- scatterplot.wrapper(
            data = baf.logr.per.chr,
            x = baf.logr.per.chr$Position,
            y = logr.name,
            xaxis.label = NULL,
            xlab.position = TRUE,
            xgrid.position = xgrid.position,
            ylab.position = ylab.position.logr,
            xlab = NULL,
            ylab = ylab.logr,
            ylim = ylim.logr,
            abline.pos = abline.pos.logr
            );
        generate.multipanelplot(
            plot.file.path = plot.file.path,
            plot.title = main.title
            );
        }
    }

### Full Dataset - Plot BAF and logR based on the x-axis style
if (!isTRUE(per.chrom.mode)) {
    print('Generating LogR/BAF plot on full dataset!');
    main.title <- sample.id;
    plot.file.path <- file.path(
        args$output_dir,
        paste0(sample.id, '_LogR-BAF-', variant.position.scale, '.png')
        );

    if ('index' == variant.position.scale) {
        chr.grid.position <- naturalsort(tapply(1:length(baf.logr$Chromosome), baf.logr$Chromosome, max));
        chr.grid.position <- c(0, chr.grid.position);
        plot.xlabel.position <- generate.xlabel.position(
            xgrid.position = chr.grid.position
            );

        baf.plot <- scatterplot.wrapper(
            data = baf.logr,
            x = 1:nrow(baf.logr),
            y = baf.name,
            xaxis.label = sub('chr', '', chr.unique),
            xlab.position = plot.xlabel.position,
            xgrid.position = chr.grid.position,
            ylab.position = ylab.position.baf,
            xlab = 'Chromosome',
            ylab = ylab.baf,
            ylim = ylim.baf,
            abline.pos = abline.pos.baf
            );

        logr.plot <- scatterplot.wrapper(
            data = baf.logr,
            x = 1:nrow(baf.logr),
            y = logr.name,
            xaxis.label = NULL,
            xlab.position = plot.xlabel.position,
            xgrid.position = chr.grid.position,
            ylab.position = ylab.position.logr,
            xlab = NULL,
            ylab = ylab.logr,
            ylim = ylim.logr,
            abline.pos = abline.pos.logr
            );
        } else if ('genome-position' == variant.position.scale) {
            ### Convert variant positions to kbp
            baf.logr$Position <- basepairs.to.kilobasepairs(x.basepairs = baf.logr$Position);
            ### Spatially arrange variants based on their position in the genome
            data.position.modified <- modify.positions(baf.logr, chr.length);
            baf.logr <- data.position.modified$df.mod.pos;
            chr.grid.position <- c(0, data.position.modified$grid.pos);
            plot.xlabel.position <- generate.xlabel.position(
                xgrid.position = chr.grid.position
                );

            baf.plot <- scatterplot.wrapper(
                data = baf.logr,
                x = baf.logr$Position.mod,
                y = baf.name,
                xaxis.label = sub('chr', '', chr.unique),
                xlab.position = plot.xlabel.position,
                xgrid.position = chr.grid.position,
                ylab.position = ylab.position.baf,
                xlab = 'Chromosome',
                ylab = ylab.baf,
                ylim = ylim.baf,
                abline.pos = abline.pos.baf
                );

            logr.plot <- scatterplot.wrapper(
                data = baf.logr,
                x = baf.logr$Position.mod,
                y = logr.name,
                xaxis.label = NULL,
                xlab.position = plot.xlabel.position,
                xgrid.position = chr.grid.position,
                ylab.position = ylab.position.logr,
                xlab = NULL,
                ylab = ylab.logr,
                ylim = ylim.logr,
                abline.pos = abline.pos.logr
                );
            }
    generate.multipanelplot(
        plot.file.path = plot.file.path,
        plot.title = main.title
        );
    }
