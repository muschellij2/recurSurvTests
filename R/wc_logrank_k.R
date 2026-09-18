#' K-group Luo-Huang weighted-risk-set recurrent gap-time rank test
#'
#' Extends the two-group [wc_logrank()] score-vector construction to `K >= 2`
#' groups. It uses subject-level influence vectors and their empirical covariance
#' so repeated gaps remain clustered. For two groups it is an omnibus version of
#' [wc_logrank()].
#' @param data,group,id,gap,status,episode,rho,tau As for [wc_logrank()].
#' @return An omnibus chi-square test, score vector, covariance, and event-time details.
#' @references Luo X, Huang CY (2011). Analysis of recurrent gap time data using
#' the weighted risk-set method. *Statistics in Medicine*, 30, 301-311.
#' @export
wc_logrank_k <- function(data, group, id="id", gap="gap", status="status", episode=NULL, rho=0, tau=Inf) {
  d <- .validate_gap_data(data,id,gap,status,episode); if (!group %in% names(d)) stop("Missing group column: ",group)
  recs <- .subject_records(d,id,gap,status,group=group); lev <- unique(vapply(recs,function(x) as.character(x$group),"")); if(length(lev)<2) stop("at least two groups required")
  tt <- .event_grid(recs,tau); if(!length(tt)) stop("No completed recurrent gaps within tau")
  cc <- .wrs_components(recs,tt); R <- colSums(cc$R_i); dN <- colSums(cc$dN_i); dl <- dN/R
  sl <- c(1, utils::head(cumprod(pmax(0,1-dl)),-1)); w <- sl^rho
  Z <- sapply(lev[-1], function(g) as.numeric(vapply(recs,function(x)as.character(x$group),"")==g)); if(is.null(dim(Z))) Z <- matrix(Z,ncol=1)
  psi <- matrix(0,nrow(Z),ncol(Z)); score <- numeric(ncol(Z))
  for(j in seq_along(tt)) { e <- colSums(sweep(Z,1,cc$R_i[,j],"*"))/R[j]; ze <- sweep(Z,2,e,"-"); score <- score + w[j]*colSums(sweep(ze,1,cc$dN_i[,j],"*")); dm <- cc$dN_i[,j]-cc$R_i[,j]*dl[j]; psi <- psi + w[j]*sweep(ze,1,dm,"*") }
  G <- score/nrow(Z); V <- crossprod(psi)/nrow(Z)^2; pin <- .pinv_rank(V); stat <- if(pin$rank) as.numeric(t(G)%*%pin$inv%*%G) else NA_real_
  list(method="K-group Luo-Huang WRS G_rho* omnibus test (score-vector extension)", rho=rho, groups=lev, score=stats::setNames(G,lev[-1]), cov=V, chisq=stat, df=pin$rank, p.value=if(is.na(stat)) NA_real_ else stats::pchisq(stat,pin$rank,lower.tail=FALSE), detail=data.frame(time=tt,risk=R,dN=dN,weight=w))
}
