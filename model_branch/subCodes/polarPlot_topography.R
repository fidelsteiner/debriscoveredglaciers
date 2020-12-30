polarPlot_topography <- function (mydata, pollutant = "nox", x = "ws", wd = "wd", type = "default", 
          statistic = "mean", resolution = "fine", limits = NA, exclude.missing = TRUE, 
          uncertainty = FALSE, percentile = NA, cols = "default", 
          weights = c(0.25, 0.5, 0.75), min.bin = 1, mis.col = "grey", 
          alpha = 1, upper = NA, angle.scale = 315, units = x, force.positive = TRUE, 
          k = 100, normalise = FALSE, key.header = "", key.footer = pollutant, 
          key.position = "right", key = TRUE, auto.text = TRUE, ws_spread = 15, 
          wd_spread = 4, kernel = "gaussian", tau = 0.5, ...) 
{
  library('openair')
  z <- . <- NULL
  if (statistic == "percentile" & is.na(percentile[1] & statistic != 
                                        "cpf")) {
    warning("percentile value missing, using 50")
    percentile <- 50
  }
  statistic <- gsub("\\.| ", "_", statistic)
  correlation_stats <- c("r", "slope", "intercept", "robust_slope", 
                         "robust_intercept", "quantile_slope", "quantile_intercept")
  if (statistic %in% correlation_stats && length(pollutant) != 
      2) {
    stop("Correlation statistic requires two pollutants.")
  }
  nam.x <- x
  nam.wd <- wd
  if (length(type) > 2) 
    stop("Maximum number of types is 2.")
  if (uncertainty) 
    type <- "default"
  if (uncertainty & length(pollutant) > 1) {
    stop("Can only have one pollutant when uncertainty = TRUE")
  }
  if (!statistic %in% c("mean", "median", "frequency", "max", 
                        "stdev", "weighted_mean", "percentile", "cpf", correlation_stats)) {
    stop(paste0("statistic '", statistic, "' not recognised."), 
         call. = FALSE)
  }
  if (length(weights) != 3) 
    stop("weights should be of length 3.")
  if (missing(key.header)) 
    key.header <- statistic
  if (key.header == "weighted_mean") 
    key.header <- "weighted\nmean"
  if (key.header == "percentile") {
    key.header <- c(paste(percentile, "th", sep = ""), "percentile")
  }
  if ("cpf" %in% key.header) 
    key.header <- c("CPF", "probability")
  if (length(cols) == 1 && cols == "greyscale") {
    trellis.par.set(list(strip.background = list(col = "white")))
  }
  current.strip <- trellis.par.get("strip.background")
  current.font <- trellis.par.get("fontsize")
  on.exit(trellis.par.set(fontsize = current.font))
  extra.args <- list(...)
  extra.args$xlab <- if ("xlab" %in% names(extra.args)) {
    quickText(extra.args$xlab, auto.text)
  }
  else {
    quickText("", auto.text)
  }
  extra.args$ylab <- if ("ylab" %in% names(extra.args)) {
    quickText(extra.args$ylab, auto.text)
  }
  else {
    quickText("", auto.text)
  }
  extra.args$main <- if ("main" %in% names(extra.args)) {
    quickText(extra.args$main, auto.text)
  }
  else {
    quickText("", auto.text)
  }
  if ("fontsize" %in% names(extra.args)) {
    trellis.par.set(fontsize = list(text = extra.args$fontsize))
  }
  if (!"layout" %in% names(extra.args)) {
    extra.args$layout <- NULL
  }
  vars <- c(wd, x, pollutant)
  if (any(type %in% dateTypes)) 
    vars <- c(vars, "date")
  mydata <- checkPrep(mydata, vars, type, remove.calm = FALSE)
  mydata <- na.omit(mydata)
  min.scale <- min(mydata[[x]], na.rm = TRUE)
  if (length(which(mydata[pollutant] < 0))/nrow(mydata) > 
      0.1 && force.positive) {
    warning(">10% negative data detected, set force.positive = FALSE?")
  }
  mydata[[x]] <- mydata[[x]] - min(mydata[[x]], na.rm = TRUE)
  if (length(pollutant) > 1 && !statistic %in% correlation_stats) {
    if (length(type) > 1) {
      warning(paste("Only type = '", type[1], "' will be used", 
                    sep = ""))
      type <- type[1]
    }
    mydata <- gather(mydata, key = variable, value = value, 
                     UQS(syms(pollutant)), factor_key = TRUE)
    pollutant <- "value"
    if (type == "default") {
      type <- "variable"
    }
    else {
      type <- c(type, "variable")
    }
  }
  mydata <- cutData(mydata, type, ...)
  max.ws <- max(mydata[[x]], na.rm = TRUE)
  min.ws <- min(mydata[[x]], na.rm = TRUE)
  clip <- TRUE
  if (missing(upper)) {
    upper <- max.ws
    clip <- FALSE
  }
  if (resolution == "normal") 
    int <- 101
  if (resolution == "fine") 
    int <- 201
  if (resolution == "ultra.fine") 
    int <- 401
  if (all(mydata[[wd]]%%10 == 0, na.rm = TRUE)) {
    wd.int <- 10
  }
  else {
    wd.int <- 5
  }
  ws.seq <- seq(min.ws, max.ws, length = 30)
  wd.seq <- seq(from = wd.int, to = 360, by = wd.int)
  ws.wd <- expand.grid(x = ws.seq, wd = wd.seq)
  u <- with(ws.wd, x * sin(pi * wd/180))
  v <- with(ws.wd, x * cos(pi * wd/180))
  input.data <- expand.grid(u = seq(-upper, upper, length = int), 
                            v = seq(-upper, upper, length = int))
  if (statistic == "cpf") {
    if (length(percentile) > 1) {
      statistic <- "cpfi"
      if (length(percentile) == 3) {
        Mean <- mean(mydata[[pollutant]], na.rm = TRUE)
        if (percentile[3] < 0) {
          Pval <- percentile[1:2]
        }
        else {
          Pval <- quantile(subset(mydata[[pollutant]], 
                                  mydata[[pollutant]] >= Mean * percentile[3]), 
                           probs = percentile[1:2]/100, na.rm = TRUE)
        }
      }
      else {
        Pval <- quantile(mydata[[pollutant]], probs = percentile/100, 
                         na.rm = TRUE)
      }
      sub <- paste("CPF (", format(Pval[1], digits = 2), 
                   " to ", format(Pval[2], digits = 2), ")", sep = "")
    }
    else {
      Pval <- quantile(mydata[[pollutant]], probs = percentile/100, 
                       na.rm = TRUE)
      sub <- paste("CPF at the ", percentile, "th percentile (=", 
                   format(Pval, digits = 2), ")", sep = "")
    }
  }
  else {
    sub <- NULL
  }
  prepare.grid <- function(mydata) {
    wd <- cut(wd.int * ceiling(mydata[[wd]]/wd.int - 0.5), 
              breaks = seq(0, 360, wd.int), include.lowest = TRUE)
    x <- cut(mydata[[x]], breaks = seq(0, max.ws, length = 31), 
             include.lowest = TRUE)
    if (!statistic %in% correlation_stats) {
      binned <- switch(statistic, frequency = tapply(mydata[[pollutant]], 
                                                     list(wd, x), function(x) length(na.omit(x))), 
                       mean = tapply(mydata[[pollutant]], list(wd, 
                                                               x), function(x) mean(x, na.rm = TRUE)), median = tapply(mydata[[pollutant]], 
                                                                                                                       list(wd, x), function(x) median(x, na.rm = TRUE)), 
                       max = tapply(mydata[[pollutant]], list(wd, x), 
                                    function(x) max(x, na.rm = TRUE)), stdev = tapply(mydata[[pollutant]], 
                                                                                      list(wd, x), function(x) sd(x, na.rm = TRUE)), 
                       cpf = tapply(mydata[[pollutant]], list(wd, x), 
                                    function(x) (length(which(x > Pval))/length(x))), 
                       cpfi = tapply(mydata[[pollutant]], list(wd, 
                                                               x), function(x) (length(which(x > Pval[1] & 
                                                                                               x <= Pval[2]))/length(x))), weighted_mean = tapply(mydata[[pollutant]], 
                                                                                                                                                  list(wd, x), function(x) (mean(x) * length(x)/nrow(mydata))), 
                       percentile = tapply(mydata[[pollutant]], list(wd, 
                                                                     x), function(x) quantile(x, probs = percentile/100, 
                                                                                              na.rm = TRUE)))
      binned <- as.vector(t(binned))
    }
    else {
      binned <- rowwise(ws.wd) %>% do(calculate_weighted_statistics(., 
                                                                    mydata, statistic = statistic, x = nam.x, y = nam.wd, 
                                                                    pol_1 = pollutant[1], pol_2 = pollutant[2], 
                                                                    ws_spread = ws_spread, wd_spread = wd_spread, 
                                                                    kernel, tau = tau))
      binned <- binned$stat_weighted
      binned <- ifelse(binned == Inf, NA, binned)
    }
    bin.len <- tapply(mydata[[pollutant[1]]], list(x, wd), 
                      length)
    binned.len <- as.vector(bin.len)
    W <- rep(1, times = length(binned))
    ids <- which(binned.len == 1)
    W[ids] <- W[ids] * weights[1]
    ids <- which(binned.len == 2)
    W[ids] <- W[ids] * weights[2]
    ids <- which(binned.len == 3)
    W[ids] <- W[ids] * weights[3]
    ids <- which(binned.len < min.bin)
    binned[ids] <- NA
    binned.len[ids] <- NA
    if (force.positive) 
      n <- 0.5
    else n <- 1
    if (!uncertainty) {
      Mgam <- try(gam(binned^n ~ s(u, v, k = k), weights = W), 
                  TRUE)
      if (!inherits(Mgam, "try-error")) {
        pred <- predict.gam(Mgam, input.data)
        pred <- pred^(1/n)
        pred <- as.vector(pred)
        results <- data.frame(u = input.data$u, v = input.data$v, 
                              z = pred)
      }
      else {
        results <- data.frame(u = u, v = v, z = binned)
        exclude.missing <- FALSE
        warning(call. = FALSE, paste("Not enough data to fit surface.\nTry reducing the value of the smoothing parameter, k to less than ", 
                                     k, ".", sep = ""))
      }
    }
    else {
      Mgam <- gam(binned^n ~ s(u, v, k = k), weights = binned.len)
      pred <- predict.gam(Mgam, input.data, se.fit = TRUE)
      uncer <- 2 * as.vector(pred[[2]])
      pred <- as.vector(pred[[1]])^(1/n)
      Mgam <- gam(binned^n ~ s(u, v, k = k))
      pred <- predict.gam(Mgam, input.data)
      pred <- as.vector(pred)
      Lower <- (pred - uncer)^(1/n)
      Upper <- (pred + uncer)^(1/n)
      pred <- pred^(1/n)
      n <- length(pred)
      results <- data.frame(u = rep(input.data$u, 3), 
                            v = rep(input.data$v, 3), z = c(pred, Lower, 
                                                            Upper), default = rep(c("prediction", "lower uncertainty", 
                                                                                    "upper uncertainty"), each = n))
    }
    exclude <- function(results) {
      x <- seq(-upper, upper, length = int)
      y <- x
      res <- int
      wsp <- rep(x, res)
      wdp <- rep(y, rep(res, res))
      all.data <- na.omit(data.frame(u, v, binned.len))
      ind <- with(all.data, exclude.too.far(wsp, wdp, 
                                            u, v, dist = 0.05))
      results$z[ind] <- NA
      results
    }
    if (exclude.missing) 
      results <- exclude(results)
    results
  }
  if (!missing(min.bin)) {
    tmp <- min.bin
    min.bin <- 0
    res1 <- group_by(mydata, UQS(syms(type))) %>% do(prepare.grid(.))
    min.bin <- tmp
    res <- group_by(mydata, UQS(syms(type))) %>% do(prepare.grid(.))
    res$miss <- res1$z
  }
  else {
    res <- group_by(mydata, UQS(syms(type))) %>% do(prepare.grid(.))
  }
  if (any(res$z > 1, na.rm = TRUE) & statistic %in% c("cpf", 
                                                      "cpfi")) {
    id <- which(res$z > 1)
    res$z[id] <- 1
  }
  if (clip) 
    res$z[(res$u^2 + res$v^2)^0.5 > upper] <- NA
  strip.dat <- strip.fun(res, type, auto.text)
  strip <- strip.dat[[1]]
  strip.left <- strip.dat[[2]]
  pol.name <- strip.dat[[3]]
  if (uncertainty) 
    strip <- TRUE
  if (normalise) {
    res <- mutate(res, z = z/mean(z, na.rm = TRUE))
    if (missing(key.footer)) 
      key.footer <- "normalised\nlevel"
  }
  if (statistic == "r") {
    if (missing(key.footer)) {
      key.footer <- paste0("corr(", pollutant[1], ", ", 
                           pollutant[2], ")")
    }
    id <- which(res$z > 1)
    if (length(id) > 0) 
      res$z[id] <- 1
    id <- which(res$z < -1)
    if (length(id) > 0) 
      res$z[id] <- -1
  }
  if (statistic == "r") 
    key.header <- expression(italic("r"))
  if (statistic == "robust_slope") 
    key.header <- "robust\nslope"
  if (statistic == "robust_intercept") 
    key.header <- "robust\nintercept"
  if (statistic == "quantile_slope") {
    key.header <- paste0("quantile slope\n(tau: ", tau, 
                         ")")
  }
  if (statistic == "quantile_intercept") {
    key.header <- paste0("quantile intercept\n(tau: ", tau, 
                         ")")
  }
  nlev <- 200
  if (missing(limits)) {
    breaks <- seq(min(res$z, na.rm = TRUE), max(res$z, na.rm = TRUE), 
                  length.out = nlev)
    labs <- pretty(breaks, 7)
    labs <- labs[labs >= min(breaks) & labs <= max(breaks)]
    at <- labs
  }
  else {
    breaks <- seq(min(limits), max(limits), length.out = nlev)
    labs <- pretty(breaks, 7)
    labs <- labs[labs >= min(breaks) & labs <= max(breaks)]
    at <- labs
    if (max(limits) < max(res[["z"]], na.rm = TRUE)) {
      id <- which(res[["z"]] > max(limits))
      res[["z"]][id] <- max(limits)
      labs[length(labs)] <- paste(">", labs[length(labs)])
    }
    if (min(limits) > min(res[["z"]], na.rm = TRUE)) {
      id <- which(res[["z"]] < min(limits))
      res[["z"]][id] <- min(limits)
      labs[1] <- paste("<", labs[1])
    }
  }
  nlev2 <- length(breaks)
  col <- openColours(cols, (nlev2 - 1))
  col.scale <- breaks
  if (uncertainty & is.null(extra.args$layout)) 
    extra.args$layout <- c(3, 1)
  legend <- list(col = col, at = col.scale, labels = list(labels = labs, 
                                                          at = at), space = key.position, auto.text = auto.text, 
                 footer = key.footer, header = key.header, height = 1, 
                 width = 1.5, fit = "all")
  legend <- makeOpenKeyLegend(key, legend, "polarPlot")
  intervals <- pretty(c(mydata[[x]], upper))
  labels <- pretty(c(mydata[[x]], upper) + min.scale)
  intervals <- intervals + (min(labels) - min.scale)
  if (min.scale != 0) {
    labels <- labels[-1]
    intervals <- intervals[-1]
  }
  temp <- paste(type, collapse = "+")
  myform <- formula(paste("z ~ u * v | ", temp, sep = ""))
  Args <- list(x = myform, res, axes = FALSE, as.table = TRUE, 
               strip = strip, strip.left = strip.left, col.regions = col, 
               region = TRUE, aspect = 1, sub = sub, par.strip.text = list(cex = 0.8), 
               scales = list(draw = FALSE), xlim = c(-upper * 1.025, 
                                                     upper * 1.025), ylim = c(-upper * 1.025, upper * 
                                                                                1.025), colorkey = FALSE, legend = legend, panel = function(x, 
                                                                                                                                            y, z, subscripts, ...) {
                                                                                  if (min.bin > 1) {
                                                                                    panel.levelplot(x, y, res$miss, subscripts, 
                                                                                                    col.regions = mis.col, labels = FALSE)
                                                                                  }
                                                                                  panel.levelplot(x, y, z, subscripts, at = col.scale, 
                                                                                                  pretty = TRUE, col.regions = col, labels = FALSE, 
                                                                                                  alpha.regions = alpha)
                                                                                  angles <- seq(0, 2 * pi, length = 360)
                                                                                  sapply(intervals, function(x) llines(x * sin(angles), 
                                                                                                                       x * cos(angles), col = "grey", lty = 5))
                                                                                  ltext(1.07 * intervals * sin(pi * angle.scale/180), 
                                                                                        1.07 * intervals * cos(pi * angle.scale/180), 
                                                                                        sapply(paste(labels, c("", "", units, rep("", 
                                                                                                                                  7))), function(x) quickText(x, auto.text)), 
                                                                                        cex = 0.7, pos = 4)
                                                                                  lsegments(-upper, 0, upper, 0)
                                                                                  lsegments(0, -upper, 0, upper)
                                                                                  ltext(upper * -1 * 0.95, 0.07 * upper, "W", cex = 0.7)
                                                                                  ltext(0.07 * upper, upper * -1 * 0.95, "S", cex = 0.7)
                                                                                  ltext(0.07 * upper, upper * 0.95, "N", cex = 0.7)
                                                                                  ltext(upper * 0.95, 0.07 * upper, "E", cex = 0.7)
                                                                                  if (grepl("slope|intercept", statistic) & length(pollutant == 
                                                                                                                                   2)) {
                                                                                    label_formula <- quickText(paste0("Formula:\n", 
                                                                                                                      pollutant[1], " ~ ", pollutant[2]))
                                                                                    ltext(upper * 0.8, 0.8 * upper, label_formula, 
                                                                                          cex = 0.7)
                                                                                  }
                                                                                })
  Args <- listUpdate(Args, extra.args)
  plt <- do.call(levelplot, Args)
  if (length(type) == 1) {
    plot(plt)
  }
  else {
    plot(useOuterStrips(plt, strip = strip, strip.left = strip.left))
  }
  newdata <- res
  output <- list(plot = plt, data = newdata, call = match.call())
  class(output) <- "openair"
  invisible(output)
}