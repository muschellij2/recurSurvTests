# Validate the K-group WRS score-vector extension; never run during checks.
# RUN_WC_LOGRANK_K_VALIDATION=true Rscript data-raw/wc-logrank-k-validation-simulation.R
if (Sys.getenv("RUN_WC_LOGRANK_K_VALIDATION") != "true") {
  message("Set RUN_WC_LOGRANK_K_VALIDATION=true to run.")
} else {
  devtools::load_all(".", quiet = TRUE)
  n_total <- as.integer(Sys.getenv("WC_LOGRANK_K_N_SIM", "5000"))
  task_id <- as.integer(Sys.getenv("WC_LOGRANK_K_TASK_ID", "1")); n_tasks <- as.integer(Sys.getenv("WC_LOGRANK_K_N_TASKS", "1"))
  n_sim <- ceiling(n_total / n_tasks)
  set.seed(as.integer(Sys.getenv("WC_LOGRANK_K_SEED", "20260923")) + task_id - 1L)
  # PSG-shaped: unequal arms, lognormal shared propensity, burst/recovery gaps,
  # 360--480 minute terminal censoring, and roughly 200 bouts per visit.
  one_subject <- function(id, arm, multiplier) {
    total <- 0; state <- stats::rbinom(1,1,.5); frailty <- stats::rlnorm(1,0,.55)
    end <- floor(stats::runif(1,360,480)/.5)*.5; g <- numeric()
    repeat { y <- ceiling(stats::rlnorm(1, log(2.4*frailty*multiplier*(if(state) .55 else 1.8)), .4)/.5)*.5
      if(total+y >= end) { g <- c(g,end-total); break }; g <- c(g,y); total <- total+y
      if(stats::runif(1) > .85) state <- 1-state }
    data.frame(id=id, arm=arm, episode=seq_along(g), gap=g, status=c(rep(1L,length(g)-1),0L))
  }
  simulate <- function(mult=c(A=1,B=1,C=1), n=c(A=30,B=45,C=60)) {
    a <- rep(names(n),n); do.call(rbind,lapply(seq_along(a),function(i) one_subject(i,a[i],mult[a[i]])))
  }
  scenarios <- list(null=c(A=1,B=1,C=1), ordered=c(A=1,B=1.2,C=.8))
  ans <- do.call(rbind,lapply(names(scenarios),function(s) {
    pv <- replicate(n_sim, wc_logrank_k(simulate(scenarios[[s]]),group="arm",episode="episode")$p.value)
    rate <- mean(pv < .05); data.frame(scenario=s,rejection_rate=rate,mc_se=sqrt(rate*(1-rate)/n_sim),n_sim=n_sim)
  }))
  suffix <- if (n_tasks > 1L) paste0("-task", task_id) else ""
  utils::write.csv(ans,paste0("data-raw/wc-logrank-k-validation-results",suffix,".csv"),row.names=FALSE)
  saveRDS(ans,paste0("data-raw/wc-logrank-k-validation-results",suffix,".rds")); print(ans)
}
