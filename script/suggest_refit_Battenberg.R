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
    '-f',
    '--default-refit-suggestions',
    type = 'character',
    help = 'Path to Battenberg  default `*refit-suggestion.txt` file for the tumor sample',
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
    '-r',
    '--min-rho',
    type = 'double',
    help = 'Min tumor purity (rho) threshold for a CNA segment',
    default = 0.1,
    required = TRUE
    );
parser$add_argument(
    '-m',
    '--min-psi',
    type = 'double',
    help = 'Min tumor ploidy (psi) threshold for a CNA segment',
    default = 1.5,
    required = TRUE
    );
parser$add_argument(
    '-n',
    '--max-psi',
    type = 'double',
    help = 'Max tumor ploidy (psi) threshold for a CNA segment',
    default = 2.7,
    required = TRUE
    );
parser$add_argument(
    '-g',
    '--platform-gamma',
    type = 'double',
    help = 'Platform specific gamma value (0.55 for SNP6, 1 for NGS)',
    default = 1,
    required = TRUE
    );
parser$add_argument(
    '-o',
    '--output-dir',
    type = 'character',
    help = 'Output dir path to publish refit suggestions for Battenberg output',
    required = TRUE
    );
parser$add_argument(
    '-u',
    '--util-script',
    type = 'character',
    help = 'Path to refit utils R script',
    required = TRUE
    );

### Parse arguments
args <- parser$parse_args();

tumor.id <- args$tumor_id;
subclones.file <- args$subclones_file;
default.refit.suggestions <- args$default_refit_suggestions;
min.rho <- args$min_rho;
min.psi <- args$min_psi;
max.psi <- args$max_psi;
platform.gamma <- args$platform_gamma;
output.dir <- args$output_dir;
util.script <- args$util_script;

### load Refit Utilities
source(util.script);

subclones <- read.delim(subclones.file);
default.suggest.refit.file <- default.refit.suggestions;
default.suggest.refit <- read.delim(default.suggest.refit.file);

### Run custom refit suggestions
print(
    paste(
        'Running custom refit suggestions with:',
        paste('min.rho =', min.rho),
        paste('min.psi =', min.psi),
        paste('max.psi =', max.psi),
        sep = ' '
        )
    );
custom.suggest.refit.file <- file.path(output.dir, paste0(tumor.id, '_custom_refit_suggestion.txt'));

custom.suggest.refit <- custom.refit.suggestion(
    subclones,
    min.rho,
    min.psi,
    max.psi
    );

if (is.null(validate.suggestion(custom.suggest.refit))) {
    print('No custom refit suggestions were made!')
    } else {
        write.table(
            x = custom.suggest.refit,
            file = custom.suggest.refit.file,
            sep = '\t',
            quote = FALSE,
            row.names = FALSE
            );
        print(paste('User defined refit suggestions were generated', custom.suggest.refit.file));
        }

### Validate suggest refit files
print(paste('Validating default refit suggestions file,', default.suggest.refit.file));

if (is.null(validate.suggestion(default.suggest.refit))) {
    print(
        paste(
            'No refit suggestions in default file,',
            default.suggest.refit.file
            )
        );
    if (file.exists(custom.suggest.refit.file)) {
        print(
            paste(
                'Use custom refit suggestions in,',
                custom.suggest.refit.file
                )
            );
        }
    } else {
        print('Default refit suggestions validated!');
    }
