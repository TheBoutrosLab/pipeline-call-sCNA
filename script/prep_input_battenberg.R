### Required packages
required.packages <- c(
    'argparse'
    );
### Check and install missing packages
missing.packages <- required.packages[!(required.packages %in% installed.packages())];
if (length(missing.packages)) {
    install.packages(missing.packages);
    }

### Load required packages
lapply(
    X = required.packages,
    FUN = library,
    character.only = TRUE
    );

### Obtain command line arguments
parser <- ArgumentParser();
parser$add_argument(
    '-s',
    '--subclones-file',
    type = 'character',
    help = 'Path to Battenberg `*subclones.txt` file for the tumor sample',
    required = TRUE
    );

parser$add_argument(
    '-t',
    '--tumor-id',
    type = 'character',
    help = 'Tumor sample ID',
    required = TRUE
    );

parser$add_argument(
    '-o',
    '--output-dir',
    type = 'character',
    help = 'Output dir path to store Battenberg input for CNA signatures workflow',
    required = TRUE
    );

### Parse arguments
args <- parser$parse_args();

tumor.id <- args$tumor_id;
subclones.file <- args$subclones_file;
output.dir <- args$output_dir;

subclones <- read.delim(subclones.file);
subclones$sample <- tumor.id;

# Reorder and select desired columns
clonal <- subclones[, c('sample', 'chr', 'startpos', 'endpos', 'nMaj1_A', 'nMin1_A', 'nMaj2_A', 'nMin2_A')];

write.table(
    clonal,
    file = file.path(output.dir, paste0(tumor.id, '_seg_input.txt')),
    sep = '\t',
    quote = FALSE,
    row.names = FALSE
    );
