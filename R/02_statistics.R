
valve_p_display <- function(p, method = "") {
  if (length(method) == 1L) method <- rep(method, length(p))
  ans <- ifelse(is.na(p), "NE", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
  mark <- !is.na(p) & !is.na(method) & grepl("Firth", method, fixed = TRUE)
  ans[mark] <- paste0(ans[mark], "\u2020")
  ans
}

valve_capture <- function(expr) {
  messages <- character()
  value <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = function(e) {
    messages <<- c(messages, conditionMessage(e)); NULL
  })
  list(value = value, warnings = paste(unique(messages), collapse = "; "))
}

valve_wilson <- function(events, n, conf.level = 0.95) {
  if (length(events) != 1L || length(n) != 1L || is.na(events) || is.na(n) ||
      n < 0 || events < 0 || events > n) stop("Invalid Wilson events/denominator.")
  events <- unname(events); n <- unname(n)
  if (n == 0) return(c(Proportion = NA_real_, Lower = NA_real_, Upper = NA_real_))
  z <- qnorm((1 + conf.level)/2); p <- events/n; den <- 1 + z^2/n
  centre <- (p + z^2/(2*n))/den
  half <- z*sqrt(p*(1-p)/n + z^2/(4*n^2))/den
  c(Proportion = p, Lower = max(0, centre-half), Upper = min(1, centre+half))
}

valve_fisher <- function(m, config = list(), context = "") {
  # Rows are groups/subgroup categories; columns are outcome categories.
  # Empty specified groups remain NE; empty outcome columns may be removed.
  m <- as.matrix(m)
  method <- if (nrow(m) == 2L && ncol(m) == 2L) {
    "Two-sided Fisher exact test (2x2)"
  } else sprintf("Global two-sided Fisher-Freeman-Halton exact test (%dx%d)", nrow(m), ncol(m))
  out <- list(P = NA_real_, Method = method, Status = "OK", Warnings = "")
  if (anyNA(m) || any(!is.finite(m)) || any(m < 0) || any(m != floor(m))) {
    out$Status <- "NE: invalid contingency-table counts"; return(out)
  }
  if (nrow(m) < 2L || any(rowSums(m) == 0)) {
    out$Status <- "NE: empty group/subgroup denominator"; return(out)
  }
  m <- m[, colSums(m) > 0, drop = FALSE]
  if (ncol(m) <= 1L) { out$P <- 1; return(out) }
  workspace <- if (is.null(config$fisher_workspace)) 2e7 else config$fisher_workspace
  test <- valve_capture(stats::fisher.test(m, alternative = "two.sided",
                         workspace = workspace, simulate.p.value = FALSE))
  out$Warnings <- test$warnings
  if (is.null(test$value)) {
    out$Status <- "NE: exact Fisher computation failed; no simulated substitution"
  } else out$P <- unname(test$value$p.value)
  out
}

valve_counts <- function(d, outcome) {
  if (!outcome %in% names(d)) stop("Outcome column missing: ", outcome)
  ans <- lapply(c("Intervention", "Control"), function(arm) {
    y <- d[[outcome]][d$randomization_group %in% arm]
    if (any(!is.na(y) & !y %in% c(0, 1))) stop("Endpoint must be 0/1/NA: ", outcome)
    data.frame(Arm = arm, N = sum(!is.na(y)), Events = sum(y == 1, na.rm = TRUE),
               Missing = sum(is.na(y)), stringsAsFactors = FALSE)
  })
  do.call(rbind, ans)
}

valve_frequency_data <- function(counts) {
  d <- do.call(rbind, lapply(seq_len(nrow(counts)), function(i) {
    data.frame(y = c(0, 1), randomization_group = counts$Arm[i],
               SUBGROUP = if ("Stratum" %in% names(counts)) counts$Stratum[i] else "All",
               freq = c(counts$N[i]-counts$Events[i], counts$Events[i]),
               stringsAsFactors = FALSE)
  }))
  d <- d[d$freq > 0, , drop = FALSE]
  d$randomization_group <- factor(d$randomization_group, levels = c("Control", "Intervention"))
  d$SUBGROUP <- factor(d$SUBGROUP)
  d
}

valve_firth_functions <- function() {
  if (requireNamespace("logistf", quietly = TRUE)) {
    return(list(fit = logistf::logistf, test = logistf::logistftest,
       control = logistf::logistf.control(maxit = 1000, collapse = FALSE),
       plcontrol = logistf::logistpl.control(maxit = 1000)))
  }
  stop("Package 'logistf' is required for zero-event blocks. Run install.packages('logistf').")
}

valve_effect <- function(counts, method = "glm") {
  out <- list(Method = method, Status = "OK", OR = NA_real_, Lower = NA_real_,
              Upper = NA_real_, P = NA_real_, Warnings = "")
  if (nrow(counts) != 2L || any(counts$N == 0)) out$Status <- "Empty arm"
  else if (all(counts$Events == 0)) out$Status <- "Both arms zero events: NE"
  if (out$Status != "OK") { out$Method <- "NE"; return(out) }
  d <- valve_frequency_data(counts)
  if (method == "glm") {
    result <- valve_capture(stats::glm(y ~ randomization_group, data = d, weights = freq,
                             family = stats::binomial(),
                             control = stats::glm.control(maxit = 100, epsilon = 1e-12)))
    fit <- result$value; out$Warnings <- result$warnings
    if (is.null(fit) || !isTRUE(fit$converged) || any(counts$N-counts$Events == 0)) {
      out$Status <- "Separated or nonconvergent glm"
    } else {
      z <- summary(fit)$coefficients["randomization_groupIntervention", ]
      out$OR <- unname(exp(z[1])); out$Lower <- unname(exp(z[1]-qnorm(.975)*z[2]))
      out$Upper <- unname(exp(z[1]+qnorm(.975)*z[2])); out$P <- unname(z[4])
    }
  } else {
    api <- valve_firth_functions()
    result <- valve_capture(api$fit(y ~ randomization_group, data = d, weights = freq,
                  firth = TRUE, pl = TRUE, plconf = 2, model = TRUE,
                  control = api$control, plcontrol = api$plcontrol))
    fit <- result$value; out$Warnings <- result$warnings
    if (is.null(fit) || any(fit$iter >= 1000) || any(fit$pl.iter >= 1000) ||
        grepl("converg|iterations.*exceeded", out$Warnings, ignore.case = TRUE)) {
      out$Status <- "Firth fit/profile failed"
    } else {
      term <- "randomization_groupIntervention"
      out$OR <- unname(exp(fit$coefficients[term]))
      out$Lower <- unname(exp(fit$ci.lower[term])); out$Upper <- unname(exp(fit$ci.upper[term]))
      out$P <- unname(fit$prob[term])
    }
  }
  if (out$Status != "OK") {
    out$Method <- "NE"; out[c("OR", "Lower", "Upper", "P")] <- rep(list(NA_real_), 4)
  }
  out
}

valve_interaction <- function(counts, method = "glm") {
  out <- list(Method = paste(method, "likelihood-ratio interaction"), Status = "OK",
              P = NA_real_, LR = NA_real_, DF = NA_real_, Warnings = "")
  if (length(unique(counts$Stratum)) < 2L) out$Status <- "Fewer than two subgroup levels: NE"
  else if (any(counts$N == 0)) out$Status <- "Empty group-by-subgroup cell"
  else if (any(vapply(split(counts, counts$Stratum), function(z) all(z$Events == 0), logical(1))))
    out$Status <- "Both arms zero events within at least one stratum: NE"
  if (out$Status != "OK") { out$Method <- "NE"; return(out) }
  d <- valve_frequency_data(counts)
  if (method == "glm") {
    a <- valve_capture(stats::glm(y ~ randomization_group*SUBGROUP, data = d, weights = freq,
          family = stats::binomial(), control = stats::glm.control(maxit = 100, epsilon = 1e-12)))
    b <- valve_capture(stats::glm(y ~ randomization_group+SUBGROUP, data = d, weights = freq,
          family = stats::binomial(), control = stats::glm.control(maxit = 100, epsilon = 1e-12)))
    out$Warnings <- paste(a$warnings, b$warnings)
    if (is.null(a$value) || is.null(b$value) || !isTRUE(a$value$converged) || !isTRUE(b$value$converged) ||
        any(counts$N-counts$Events == 0)) out$Status <- "Separated or nonconvergent glm interaction"
    else {
      out$LR <- max(0, b$value$deviance-a$value$deviance)
      out$DF <- a$value$rank-b$value$rank
      out$P <- stats::pchisq(out$LR, out$DF, lower.tail = FALSE)
    }
  } else {
    api <- valve_firth_functions()
    a <- valve_capture(api$fit(y ~ randomization_group*SUBGROUP, data = d, weights = freq,
          firth = TRUE, pl = FALSE, model = TRUE, control = api$control, plcontrol = api$plcontrol))
    out$Warnings <- a$warnings
    if (is.null(a$value) || any(a$value$iter >= 1000)) out$Status <- "Firth interaction fit failed"
    else {
      terms <- which(grepl(":", names(a$value$coefficients), fixed = TRUE))
      # Joint constrained test uses the full model's Jeffreys penalty, not a
      # subtraction of independently penalized full/reduced-model likelihoods.
      b <- valve_capture(api$test(a$value, test = terms, values = rep(0, length(terms)),
                                 control = api$control))
      out$Warnings <- paste(out$Warnings, b$warnings)
      if (is.null(b$value) || grepl("converg|iterations.*exceeded", out$Warnings, ignore.case = TRUE))
        out$Status <- "Firth interaction constrained test failed"
      else {
        out$LR <- max(0, 2*unname(b$value$loglik["full"]-b$value$loglik["null"]))
        out$DF <- length(terms); out$P <- stats::pchisq(out$LR, out$DF, lower.tail = FALSE)
      }
    }
  }
  out$Warnings <- trimws(out$Warnings)
  if (out$Status != "OK") {
    out$Method <- "NE"; out[c("P", "LR", "DF")] <- rep(list(NA_real_), 3)
  }
  out
}

valve_rbind <- function(rows) {
  if (!length(rows)) return(data.frame())
  cols <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(z) {
    for (nm in setdiff(cols, names(z))) z[[nm]] <- NA
    z[, cols, drop = FALSE]
  })
  ans <- do.call(rbind, rows); rownames(ans) <- NULL; ans
}
