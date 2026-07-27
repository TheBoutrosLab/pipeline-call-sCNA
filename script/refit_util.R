###############################################################################
# Battenberg Refit Utilities
###############################################################################
## User defined refit suggestion
### min.rho = 0.1
### min.psi = 1.5
### max.psi = 2.7
### platform.gamma = 1
custom.refit.suggestion <- function(
    subclones,
    min.rho = 0.1,
    min.psi = 1.5,
    max.psi = 2.7,
    platform.gamma = 1
    ) {
    subclones$rho_estimate <- (2 * subclones$BAF - 1) / (2 * subclones$BAF - subclones$BAF * (subclones$nMaj1_A + subclones$nMin1_A) - 1 + subclones$nMaj1_A);
    subclones$psi <- (subclones$rho_estimate * (subclones$nMaj1_A + subclones$nMin1_A) + 2 - 2 * subclones$rho_estimate) / (2^(subclones$LogR / platform.gamma));
    subclones$psi_t_estimate <- (subclones$psi - 2 * (1 - subclones$rho_estimate)) / subclones$rho_estimate;
    subclones$pos <- as.integer((subclones$endpos - subclones$startpos + 1) / 10^6);
    subclones <- subclones[order(-subclones$pos),];
    subclones$pos <- paste0(subclones$pos, 'M');

    subclones.rho.finite <- subclones[
        is.finite(subclones$rho_estimate) &
        subclones$rho_estimate > min.rho &
        subclones$psi_t_estimate > min.psi &
        subclones$psi_t_estimate < max.psi &
        subclones$chr != 'X', c('chr', 'pos', 'rho_estimate', 'psi_t_estimate')
        ];

    if (nrow(subclones.rho.finite) == 0) {
        subclones.rho.finite <- data.frame(chr = NA, pos = NA, rho_estimate = NA, psi_t_estimate = NA);
    }
    names(subclones.rho.finite)[names(subclones.rho.finite) == 'chr'] <- 'chrom';
    return(subclones.rho.finite);
    }

## Validate if there are refit suggestions at all in the first place
validate.suggestion <- function(refit.suggestion) {
    if (length(refit.suggestion$chrom) == 1 && is.na(refit.suggestion$chrom)) {
        return(NULL);
        } else {
            return('Valid');
        }
    }
