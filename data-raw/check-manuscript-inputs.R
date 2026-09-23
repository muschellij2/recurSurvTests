# Run from the repository root: Rscript data-raw/check-manuscript-inputs.R
# Uses synthetic files in a temporary directory to test report plumbing only.
here::i_am('data-raw/check-manuscript-inputs.R')
source(here::here('vignettes','manuscript-helpers.R'))
root <- here::here()
input <- file.path(root,'data-raw/zhao-2020-table-simulations')
files <- list.files(input,pattern='[.]rds$',full.names=TRUE)
# Future-data fixtures test the importer only; never use them as scientific results.
fixture <- tempfile('manuscript-fixtures-'); dir.create(fixture)
invisible(file.copy(files, fixture))
original <- manuscript_zhao(root,fixture)
missing <- original[!original$available,]
template <- readRDS(files[1])
for (j in seq_len(nrow(missing))) {
 p <- missing[j,]; x <- template
 d <- x$table1[1,]
 d$table <- p$table; d$distribution <- p$distribution; d$heterogeneity <- p$heterogeneity
 h <- switch(p$heterogeneity,low=c(.5,1.5),medium=c(.1,1.9),high=c(.01,1.99))
 d$z_low <- h[1]; d$z_high <- h[2]; d$test <- p$test
 d$n_sim <- p$paper_n_sim; d$treatment_time_multiplier <- exp(p$log_time_shift)
 d$empirical_rejection <- .025; d$monte_carlo_se <- sqrt(.025*.975/d$n_sim)
 x$table1 <- x$table2 <- x$table3 <- d[FALSE,]; x[[p$table]] <- d
 x$settings$task_id <- p$task_id; x$settings$generator_version <- 'zhao_table1_loglogistic_v2'
 saveRDS(x,file.path(fixture,sprintf('zhao-2020-table-simulations-task-%03d-of-120.rds',p$task_id)))
}
full <- manuscript_zhao(root,fixture)
stopifnot(nrow(full)==120,all(full$available),all(full$design_eligible))
stopifnot(!any(grepl('Pending|Excluded',unlist(manuscript_zhao_table(full,'table3')))))
# Legacy log-logistic output is excluded even in an otherwise complete batch.
p <- file.path(fixture,'zhao-2020-table-simulations-task-004-of-120.rds')
x <- readRDS(p); x$settings$generator_version <- NULL; saveRDS(x,p)
legacy <- manuscript_zhao(root,fixture)
stopifnot(legacy$comparison[legacy$task_id==4]=='design mismatch')
# A changed variance convention must fail, rather than pool distinct methods.
x$settings$variance_method <- 'zhao_eq6'; saveRDS(x,p)
msg <- tryCatch({manuscript_zhao(root,fixture); ''},error=conditionMessage)
stopifnot(grepl('Mixed settings',msg))
# Empty directory still yields every published row and all-pending placeholders.
empty <- tempfile('manuscript-empty-'); dir.create(empty)
d <- manuscript_zhao(root,empty)
stopifnot(nrow(d)==120,!any(d$available),nrow(manuscript_zhao_table(d,'table3'))==24)
pdf(file.path(empty,'pending.pdf')); manuscript_zhao_plot(d); dev.off()
cat('Partial/full/missing results, legacy exclusion, mixed-method rejection, and pending plot verified.\n')
