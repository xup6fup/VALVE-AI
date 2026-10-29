# Post hoc sensitivity to masking TTE ascertainment in low-risk intervention
# participants. Counts are TTE recipients, not numbers of examinations.
# Public API: valve_tte_sensitivity() calculates and returns aggregate results;
# valve_draw_tte_sensitivity() draws Ex_Data_Fig_09 on the current device.
# This module has no file I/O and performs no analysis when sourced.

valve_tte_sensitivity <- function(d, config=valve_config()) {
  # Keep simulation reproducibility local to this analysis.
  had_seed <- exists('.Random.seed', envir=.GlobalEnv, inherits=FALSE)
  previous_seed <- if (had_seed) get('.Random.seed', envir=.GlobalEnv) else NULL
  previous_rng <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(previous_rng))
    if (had_seed) assign('.Random.seed', previous_seed, envir=.GlobalEnv)
    else if (exists('.Random.seed', envir=.GlobalEnv, inherits=FALSE))
      rm('.Random.seed', envir=.GlobalEnv)
  }, add=TRUE)

  required <- c('randomization_group', 'ai_ecg_risk', 'outcome_tte_after_randomization', 'outcome_vhd')
  stopifnot(all(required %in% names(d)), !anyNA(d[required]))
  stopifnot(all(d$outcome_tte_after_randomization %in% 0:1), all(d$outcome_vhd %in% 0:1))
  stopifnot(all(d$outcome_vhd <= d$outcome_tte_after_randomization))
  I <- d$randomization_group == 'Intervention'
  C <- d$randomization_group == 'Control'
  stopifnot(all(I | C), any(I), any(C))
  ni <- sum(I); nc <- sum(C)
  ei <- sum(d$outcome_vhd[I]); ec <- sum(d$outcome_vhd[C])
  ttei <- sum(d$outcome_tte_after_randomization[I]); ttec <- sum(d$outcome_tte_after_randomization[C])
  eligible <- which(I & d$ai_ecg_risk == 'Low' & d$outcome_tte_after_randomization == 1)
  y_eligible <- d$outcome_vhd[eligible]
  B <- 10000L
  seed <- 0L
  ratios <- c(1.10, 1.05, 1.00, .95, .90)
  # Floor fractional targets so TTE use never exceeds the requested ceiling.
  target <- floor(ttec * ratios + 1e-10)
  n_mask <- ttei - target
  stopifnot(all(n_mask >= 1), all(n_mask <= length(eligible)))

  # The unadjusted binomial GLM and two-sided coefficient Wald P match the
  # original primary-outcome analysis. All possible remaining event counts
  # are positive in this analysis, so the original Firth trigger never occurs.
  fit_primary <- function(events_i) {
    z <- data.frame(
      arm = factor(c('Control','Control','Intervention','Intervention'),
                   levels = c('Control','Intervention')),
      y = c(0,1,0,1), freq = c(nc-ec,ec,ni-events_i,events_i))
    fit <- glm(y ~ arm, family = binomial(), data = z, weights = freq,
               control = glm.control(maxit = 100, epsilon = 1e-12))
    stopifnot(isTRUE(fit$converged))
    coef <- summary(fit)$coefficients['armIntervention', ]
    c(OR = unname(exp(coef[1])), P = unname(coef[4]),
      lower = unname(exp(coef[1] - qnorm(.975)*coef[2])),
      upper = unname(exp(coef[1] + qnorm(.975)*coef[2])))
  }
  # Identical event counts yield identical GLMs; cache those fits exactly.
  event_grid <- seq.int(ei-sum(y_eligible), ei)
  stopifnot(min(event_grid) > 0, ec > 0)
  lookup <- t(vapply(event_grid, fit_primary, c(OR=0,P=0,lower=0,upper=0)))
  rownames(lookup) <- event_grid
  observed <- fit_primary(ei)
  RNGkind('Mersenne-Twister', 'Inversion', 'Rejection')
  lost <- matrix(0L, nrow = B, ncol = length(ratios))
  for (j in seq_along(ratios)) {
    # Reset once before the 10,000 replicates of EACH scenario, not per replicate.
    # Reusing the seed supplies common permutations and nested masks across targets.
    set.seed(0)
    for (b in seq_len(B)) {
      permutation <- sample.int(length(eligible))
      cumulative_lost <- cumsum(y_eligible[permutation])
      lost[b,j] <- cumulative_lost[n_mask[j]]
    }
  }
  remaining <- ei-lost
  p <- matrix(lookup[match(as.vector(remaining), event_grid),'P'], nrow=B)
  or <- matrix(lookup[match(as.vector(remaining), event_grid),'OR'], nrow=B)
  success <- p < .05 & or > 1
  stopifnot(all(apply(lost,1,function(x) all(diff(x)>=0))))

  # Exact hypergeometric check of the same uniform-masking experiment.
  # This is a validation result, not a separate analysis or a clinical CI.
  exact_prob <- vapply(seq_along(ratios), function(j) {
    k <- 0:sum(y_eligible)
    pk <- dhyper(k, sum(y_eligible), length(eligible)-sum(y_eligible), n_mask[j])
    test <- lookup[match(ei-k,event_grid), ,drop=FALSE]
    sum(pk * (test[,'P'] < .05 & test[,'OR'] > 1))
  }, 0)
  q_p <- t(apply(p, 2, quantile, probs=c(.025,.25,.5,.75,.975), type=1))
  q_events <- t(apply(remaining, 2, quantile, probs=c(.025,.5,.975), type=1))
  summary <- data.frame(
    target_ratio = ratios, actual_ratio = target/ttec,
    control_tte_recipients = ttec, intervention_tte_recipients = target,
    masked_low_risk_tte_recipients = n_mask,
    intervention_n = ni, control_n = nc, control_vhd_events = ec,
    intervention_vhd_median = q_events[,2],
    intervention_vhd_p025 = q_events[,1], intervention_vhd_p975 = q_events[,3],
    p_median=q_p[,3], p_p025=q_p[,1], p_p975=q_p[,5],
    p_q25=q_p[,2], p_q75=q_p[,4],
    significant_replicates=colSums(success), total_replicates=B,
    significance_retention=colMeans(success),
    exact_hypergeometric_retention=exact_prob,
    monte_carlo_se=sqrt(colMeans(success)*(1-colMeans(success))/B))
  stopifnot(all(abs(summary$significance_retention-exact_prob) <
                  5*sqrt(exact_prob*(1-exact_prob)/B)))
  # Retain aggregate results in memory; no participant-level records are returned.
  replicates <- data.frame(replicate=rep(seq_len(B),times=length(ratios)),
    target_ratio=rep(ratios,each=B), masked_vhd=as.vector(lost),
    intervention_vhd=as.vector(remaining), OR=as.vector(or), P=as.vector(p),
    significant_benefit=as.vector(success))

  # Overall-cohort TTE utilization in the observed trial and each fixed-count target.
  # These are binomial-regression Wald CIs, not simulation percentile intervals.
  fit_tte <- function(events_i) {
    z <- data.frame(
      arm = factor(c('Control','Control','Intervention','Intervention'),
                   levels=c('Control','Intervention')),
      y=c(0,1,0,1), freq=c(nc-ttec,ttec,ni-events_i,events_i))
    fit <- glm(y ~ arm, family=binomial(), data=z, weights=freq,
               control=glm.control(maxit=100,epsilon=1e-12))
    stopifnot(isTRUE(fit$converged))
    b <- summary(fit)$coefficients['armIntervention',]
    c(OR=unname(exp(b[1])), Lower=unname(exp(b[1]-qnorm(.975)*b[2])),
      Upper=unname(exp(b[1]+qnorm(.975)*b[2])), P=unname(b[4]))
  }
  tte_counts <- c(ttei,target)
  tte_est <- t(vapply(tte_counts,fit_tte,c(OR=0,Lower=0,Upper=0,P=0)))
  tte_forest <- data.frame(
    scenario=c('Observed',sprintf('Simulated %.2f-fold',ratios)),
    target_ratio=c(NA_real_,ratios), actual_count_ratio=tte_counts/ttec,
    intervention_tte=tte_counts, intervention_n=ni,
    control_tte=ttec, control_n=nc, tte_est,
    stringsAsFactors=FALSE)

  # Equal-use scenario: every primary-outcome estimate and its model CI.
  ranked <- replicates[replicates$target_ratio == 1, , drop=FALSE]
  ranked$lower <- lookup[match(ranked$intervention_vhd,event_grid),'lower']
  ranked$upper <- lookup[match(ranked$intervention_vhd,event_grid),'upper']
  ranked <- ranked[order(ranked$OR,ranked$replicate), , drop=FALSE]
  ranked$rank <- seq_len(nrow(ranked))
  ranked$percentile_rank <- 100*ranked$rank/nrow(ranked)
  ranked$intervention_n <- ni; ranked$control_n <- nc; ranked$control_vhd <- ec
  ranked <- ranked[,c('replicate','rank','percentile_rank','target_ratio','masked_vhd',
    'intervention_vhd','intervention_n','control_vhd','control_n','OR','lower','upper','P')]
  rownames(ranked) <- NULL
  pr_rows <- ceiling(c(.05,.5)*nrow(ranked))
  annotations <- ranked[pr_rows, ,drop=FALSE]
  annotations$annotation <- c('PR 5','PR 50')
  annotations$display <- sprintf('%s: OR %.3f (95%% CI %.3f-%.3f), P = %.3f',
    annotations$annotation,annotations$OR,annotations$lower,annotations$upper,annotations$P)
  stopifnot(nrow(ranked)==B,
    all(annotations$OR == as.numeric(quantile(ranked$OR,c(.05,.5),type=1))))

  structure(list(
    settings=list(B=B, seed=seed, ratios=ratios,
      rng_kind=c('Mersenne-Twister','Inversion','Rejection'),
      primary_method='Unadjusted binomial GLM; two-sided Wald P and 95% CI',
      percentile_method='Nearest rank (empirical quantile type 1)'),
    counts=list(intervention_n=ni, control_n=nc,
      intervention_vhd=ei, control_vhd=ec,
      intervention_tte=ttei, control_tte=ttec,
      eligible_low_risk_tte=length(eligible),
      eligible_low_risk_vhd=sum(y_eligible)),
    observed=observed, summary=summary, replicates=replicates,
    ranked=ranked, annotations=annotations, forest=tte_forest
  ), class='valve_tte_sensitivity')
}

valve_draw_tte_sensitivity <- function(x, config=valve_config()) {
  stopifnot(inherits(x, 'valve_tte_sensitivity'))
  ratios <- x$settings$ratios
  B <- x$settings$B
  ni <- x$counts$intervention_n; nc <- x$counts$control_n
  ttei <- x$counts$intervention_tte; ttec <- x$counts$control_tte
  target <- x$summary$intervention_tte_recipients
  observed <- x$observed
  p <- matrix(x$replicates$P, nrow=B)
  success <- matrix(x$replicates$significant_benefit, nrow=B)
  q_p <- as.matrix(x$summary[,c('p_p025','p_q25','p_median','p_q75','p_p975')])
  ranked <- x$ranked; annotations <- x$annotations
  tte_counts <- x$forest$intervention_tte
  tte_est <- as.matrix(x$forest[,c('OR','Lower','Upper','P')])
  previous_par <- par(no.readonly=TRUE)
  on.exit({ layout(1L); par(previous_par) }, add=TRUE)

  # Paper style: black sans-serif text, blue estimates, gray table rules.
  blue <- '#3C5488'; ink <- '#000000'; muted <- '#6D6E71'
  light <- '#EEF0F2'; red <- '#E64B35'; pale_blue <- '#C8D0E1'
  positions <- seq_along(ratios)
  labels <- sprintf('%.2f',ratios)

  layout(matrix(c(1,1,2,3),nrow=2,byrow=TRUE),heights=c(2.65,4.55))
  par(oma=c(.15,.4,.15,.3),family=config$font_family,fg=ink,col.axis=ink,
      col.lab=ink,cex=.92,las=1,xaxs='i',yaxs='i',mgp=c(2.8,.6,0))
  # Panel a: TTE utilization in all randomized participants.
  par(mar=c(2.35,.6,.9,.4))
  plot(NA,xlim=c(0,1),ylim=c(.15,7.15),xaxt='n',yaxt='n',
       xlab='',ylab='',bty='n')
  rect(.005,5.55,.995,6.45,col=light,border=NA)
  text(.012,6.92,'Scenario',adj=c(0,.5),font=2,cex=.80)
  text(.285,6.92,'Intervention',font=2,cex=.80)
  text(.443,6.92,'Control',font=2,cex=.80)
  text(.840,6.92,'Odds ratio (95% CI)',font=2,cex=.78)
  text(.988,6.92,expression(italic(P)~value),adj=c(1,.5),font=2,cex=.78)
  segments(.005,6.55,.995,6.55,col=muted,lwd=.9)
  yy <- 6:1
  labels_f <- c(sprintf('Observed (%.2f-fold)',ttei/ttec),sprintf('Simulated %.2f-fold',ratios))
  text(.012,yy,labels_f,adj=c(0,.5),font=c(2,rep(1,length(ratios))),cex=.80)
  fmt_counts <- function(n,N) sprintf('%s/%s (%.1f%%)',
    format(n,big.mark=',',trim=TRUE),format(N,big.mark=',',trim=TRUE),100*n/N)
  text(.285,yy,fmt_counts(tte_counts,ni),cex=.79)
  text(.443,yy,fmt_counts(rep(ttec,length(tte_counts)),nc),cex=.79)
  f_lo <- .550; f_hi <- .755; or_min <- .5; or_max <- 2
  fx <- function(x) f_lo+(log(x)-log(or_min))/(log(or_max)-log(or_min))*(f_hi-f_lo)
  segments(fx(1),.55,fx(1),6.45,lty=2,col='#A9ADB4',lwd=.85)
  for (j in seq_along(tte_counts)) {
    segments(fx(tte_est[j,'Lower']),yy[j],fx(tte_est[j,'Upper']),yy[j],col=blue,lwd=1.3)
    points(fx(tte_est[j,'OR']),yy[j],pch=15,col=blue,cex=.65)
  }
  forest_ticks <- c(.5,.7,1,1.5,2)
  axis(1,at=fx(forest_ticks),labels=c('0.5','0.7','1','1.5','2'),
       pos=.38,tck=-.025,cex.axis=.71,col=muted,lwd=.7,padj=.15)
  mtext('Intervention / control odds ratio',side=1,line=1.10,
        at=(f_lo+f_hi)/2,cex=.65)
  text(.840,yy,sprintf('%.2f (%.2f\u2013%.2f)',tte_est[,'OR'],tte_est[,'Lower'],tte_est[,'Upper']),cex=.79)
  ptext <- ifelse(tte_est[,'P']<.001,'<0.001',sprintf('%.3f',tte_est[,'P']))
  text(.988,yy,ptext,adj=c(1,.5),cex=.79)
  mtext('a',side=3,line=.15,adj=0,font=2,cex=.95)

  # Panel b: distribution of nominal primary-outcome P values.
  par(mar=c(4.15,4.5,1.1,.7))
  plot(NA,xlim=c(.45,length(ratios)+.55),ylim=log10(c(.00020,1)),xaxt='n',yaxt='n',
       xlab='',ylab=expression('Two-sided Wald '*italic(P)*' value (log scale)'),bty='n')
  ticks <- c(.0003,.001,.003,.01,.03,.05,.1,.3,1)
  abline(h=log10(ticks),col=light,lwd=.6)
  axis(2,at=log10(ticks),labels=c('0.0003','0.001','0.003','0.01','0.03','0.05','0.10','0.30','1.00'),
       tick=FALSE,cex.axis=.78)
  for(j in positions) {
    z <- log10(p[,j]); den <- density(z,from=min(z),to=max(z),n=256)
    width <- den$y/max(den$y)*.34
    polygon(c(j-width,rev(j+width)),c(den$x,rev(den$x)),
      col=adjustcolor(blue,.28),border=blue,lwd=.8)
    q <- log10(q_p[j,])
    segments(j,q[1],j,q[5],col=blue,lwd=.9)
    segments(j-.055,q[c(1,5)],j+.055,q[c(1,5)],col=blue,lwd=.9)
    rect(j-.09,q[2],j+.09,q[4],col='white',border=blue,lwd=.9)
    segments(j-.09,q[3],j+.09,q[3],col=blue,lwd=1.4)
  }
  abline(h=log10(.05),col=red,lty=2,lwd=1.0)
  text(.53,log10(.05)+.085,expression(italic(P)=='0.05'),adj=c(0,0),col=red,cex=.74)
  abline(h=log10(observed['P']),col=ink,lty=3,lwd=1.0)
  text(3.02,log10(observed['P'])+.205,
    bquote('Observed: '*italic(P)*' = '*.(sprintf('%.6f',observed['P']))),
    adj=c(.5,0),col=ink,cex=.65)
  text(3.02,log10(observed['P'])+.075,
    sprintf('TTE = %.2f x control (%d/%d)',ttei/ttec,ttei,ttec),
    adj=c(.5,0),col=ink,cex=.65)
  text(3,log10(.91),expression('Replicates with OR > 1 and '*italic(P)<'0.05'),cex=.69,adj=c(.5,1))
  text(positions,log10(.57),sprintf('%.2f%%',colMeans(success)*100),font=2,cex=.73)
  axis(1,at=positions,labels=labels,tick=FALSE,cex.axis=.83)
  mtext(sprintf('n = %d',target),side=1,line=1.2,at=positions,cex=.65,col=ink)
  mtext('Intervention / control TTE use (fold)',side=1,line=2.85,cex=.78)
  mtext('b',side=3,line=.2,adj=-.14,font=2,cex=.95)

  # Panel c: all 10,000 equal-use estimates, sorted by OR then replicate ID.
  # Each vertical segment is one replicate's model-based Wald 95% CI.
  par(mar=c(4.15,3.8,1.1,.7))
  plot(NA,xlim=c(0,100),ylim=c(.8,2.52),xaxt='n',yaxt='n',
       xlab='',ylab='Primary-outcome odds ratio (95% CI)',bty='n')
  abline(h=c(.8,1,1.2,1.4,1.6,1.8,2,2.2),col=light,lwd=.6)
  segments(ranked$percentile_rank,ranked$lower,ranked$percentile_rank,ranked$upper,
    col=pale_blue,lwd=.18)
  points(ranked$percentile_rank,ranked$OR,pch=16,col=blue,cex=.13)
  abline(h=1,lty=2,col=muted,lwd=.85)
  axis(1,at=c(0,25,50,75,100),labels=c('0','25','50','75','100'),cex.axis=.80,
       col=muted,lwd=.7,tck=-.015)
  axis(2,at=c(.8,1,1.2,1.4,1.6,1.8,2,2.2),labels=c('0.8','1.0','1.2','1.4','1.6','1.8','2.0','2.2'),
       tick=FALSE,cex.axis=.80)
  for (j in seq_len(nrow(annotations))) {
    a <- annotations[j,]
    segments(a$percentile_rank,a$lower,a$percentile_rank,a$upper,col=blue,lwd=1.4)
    segments(a$percentile_rank-1,a$lower,a$percentile_rank+1,a$lower,col=blue,lwd=1.2)
    segments(a$percentile_rank-1,a$upper,a$percentile_rank+1,a$upper,col=blue,lwd=1.2)
    points(a$percentile_rank,a$OR,pch=21,bg=blue,col='white',cex=.9,lwd=.7)
  }
  # Gray dotted leaders distinguish annotations from the blue CI endpoints.
  segments(annotations$percentile_rank,annotations$upper+.035,
    annotations$percentile_rank,2.26,col=muted,lty=3,lwd=.6)
  segments(c(5,50),2.26,c(9,70),2.26,col=muted,lty=3,lwd=.6)
  text(9,2.47,sprintf('PR 5: OR %.3f\n95%% CI %.3f\u2013%.3f',
    annotations$OR[1],annotations$lower[1],annotations$upper[1]),
    adj=c(0,1),cex=.69,col=ink)
  text(9,2.32,bquote(italic(P)*' = '*.(sprintf('%.3f',annotations$P[1]))),
    adj=c(0,1),cex=.69,col=ink)
  text(95,2.47,sprintf('PR 50: OR %.3f\n95%% CI %.3f\u2013%.3f',
    annotations$OR[2],annotations$lower[2],annotations$upper[2]),
    adj=c(1,1),cex=.69,col=ink)
  text(95,2.32,bquote(italic(P)*' = '*.(sprintf('%.3f',annotations$P[2]))),
    adj=c(1,1),cex=.69,col=ink)
  mtext('Percentile rank of 10,000 replicates',side=1,line=2.85,cex=.78)
  mtext('c',side=3,line=.2,adj=-.12,font=2,cex=.95)
  invisible(x)
}
