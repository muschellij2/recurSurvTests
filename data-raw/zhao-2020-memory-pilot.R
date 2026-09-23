# Diagnostic only: override the log-logistic draw in an isolated environment.
# Does not modify package functions or saved simulation outputs.
# From the package root on macOS:
# /usr/bin/time -l Rscript data-raw/zhao-2020-memory-pilot.R
# On Linux, use /usr/bin/time -v instead.
# Fifteen datasets are a sizing pilot, not a worst-case memory guarantee.
e <- new.env(parent = globalenv())
for (f in list.files('R', pattern='[.]R$', full.names=TRUE)) sys.source(f, envir=e)
original_draw <- e$.zhao_draw_gap
e$.zhao_draw_gap <- function(distribution, mean_gap=NULL, time_resolution=NULL) {
 if (distribution=='loglogistic' && is.null(mean_gap)) return(exp(4)*(runif(1)^(-1)-1)^.5)
 original_draw(distribution, mean_gap, time_resolution)
}
set.seed(230923)
for (h in list(low=c(.5,1.5),medium=c(.1,1.9),high=c(.01,1.99))) {
 events <- numeric(5)
 for (i in 1:5) {
  d <- e$.zhao_simulated_gap_data(100,'loglogistic',h,180)
  events[i] <- sum(d$status)
  fit <- e$zhao_rank_test(d,group='group',episode='episode')
 }
 cat('Z range',h,'events range',range(events),'\n')
}
cat('Original-generator memory approximation (null design):\n')
for (h in list(low=c(.5,1.5),medium=c(.1,1.9),high=c(.01,1.99))) {
 M <- 200*90/(exp(-4)*pi/4)*log(h[2]/h[1])/diff(h)
 cat('Z range',h,'expected events approx',round(M),'one matrix GiB',200*M*8/1024^3,'six matrices x four workers GiB',24*200*M*8/1024^3,'\n')
}
