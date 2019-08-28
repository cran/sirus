#'
#' Fit a SIRUS classification model for a given number of rules (10 by default) or a given p0.
#' The number of trees is tuned automatically with a stopping criterion based on stability.
#' The hyperparameter p0 can be tuned using sirus.cv.
#'
#' @title Fit SIRUS.
#' @param data Input dataframe, each row is an observation vector.
#' @param y Binary response variable taking 0 and 1 values.
#' @param num.rule Number of rules in SIRUS model. Default is 10. Ignored if a p0 value is provided.
#' @param p0 Selection threshold on the frequency of appearance of a path in the forest. Default is NULL and num.rule is used to select rules.
#' @param num.rule.max Maximum number of rules in SIRUS model. Ignored if num.rule is provided.
#' @param q Number of quantiles used for node splitting in the forest construction.
#' @param num.trees.step Number of trees grown between two evaluations of the stopping criterion. Ignored if num.trees is provided.
#' @param alpha Parameter of the stopping criterion for the number of trees: stability has to reach 1 - alpha to stop the growing of the forest. Ignored if num.trees is provided.
#' @param mtry Number of variables to possibly split at each node. Default is the number of variables divided by 3.
#' @param verbose Boolean. If true, information messages are printed.
#' @param num.trees Number of trees grown in the forest. Default is NULL. If NULL (recommanded), the number of trees is automatically set using a stability based stopping criterion.
#' @param num.threads Number of threads used to grow the forest. Default is number of CPUs available.
#' @param replace Boolean. If true (default), sample with replacement.
#' @param sample.fraction Fraction of observations to sample. Default is 1 for sampling with replacement and 0.632 for sampling without replacement.
#' @param seed Random seed. Default is NULL, which generates the seed from R. Set to 0 to ignore the R seed.
#'
#' @return SIRUS model with elements
#'   \item{\code{rules}}{List of rules in SIRUS model.}
#'   \item{\code{rules.out}}{List of rule outputs. rule.out: the output mean whether the rule is satisfied or not. supp.size: the number of points inside and outside the rule.}
#'   \item{\code{proba}}{Frequency of occurence of paths in the forest.}
#'   \item{\code{paths}}{List of selected paths.}
#'   \item{\code{mean}}{Mean output over the full training data. Default model output if no rule is selected.}
#' @export
#'
#' @examples
#' ## load sirus
#' require(sirus)
#'
#' ## prepare data
#' data <- iris
#' y <- rep(0, nrow(data))
#' y[data$Species == 'setosa'] = 1
#' data$Species <- NULL
#'
#' ## fit sirus
#' sirus.m <- sirus.fit(data, y)
#'
#'
#' @encoding UTF-8
#' @useDynLib sirus
#' @importFrom Rcpp evalCpp
#' @import stats
#' @import utils
#' @importFrom Matrix Matrix
#' @import ROCR
#' @import ggplot2
sirus.fit <- function(data, y, num.rule = 10, p0 = NULL, num.rule.max = 25, q = 10,
                      num.trees.step = 1000, alpha = 0.05, mtry = NULL, num.trees = NULL,
                      num.threads = NULL, replace = TRUE, sample.fraction = ifelse(replace, 1, 0.632),
                      verbose = TRUE, seed = NULL) {

  # Check arguments
  # check data type
  data.y.check(data, y)
  # check sirus parameters
  sirus.param.check(data, num.rule.max, q, num.trees.step, alpha, mtry)
  # check num.rule
  num.rule.valid <- is.numeric(num.rule)
  if (num.rule.valid){num.rule.valid <- (round(num.rule) == num.rule) & num.rule > 0}
  if (!is.null(num.rule)){
    if (!num.rule.valid){
      stop("Invalid num.rule. Number of rules must be a positive integer or NULL.")
    }else{
      if (num.rule > 100){
        warning('Warning num.rule: SIRUS is designed to output short list of rules (typically < 100 rules).')
      }
    }
  }
  # check p0
  if (!is.null(p0)){
    p0.valid <- is.numeric(p0) & (p0 >= 0) & (p0 <= 1)
    if (!p0.valid){
      stop("Invalid p0. p0 must be a numeric value between 0 and 1 or NULL.")
    }
  }
  # check either p0 or num.rule not null
  if (is.null(p0) & is.null(num.rule)){
    stop('Invalid p0 and num.rule: Either p0 or num.rule has to be provided.')
  }
  # set default mtry
  if (is.null(mtry)){
    mtry.ratio <- 1/3
    mtry <- floor((ncol(data))*mtry.ratio)
  }

  # Data binning
  data.names <- colnames(data)
  quantiles.emp <- apply(data, 2, get.quantiles.emp, q = q)
  data.bin <- sapply(1:ncol(data), function(ind){
    get.X.bin(data[,ind], quantiles.emp[[ind]][[1]])
  })
  colnames(data.bin) <- data.names
  data.bin.y <- as.data.frame(cbind(data.bin, y))

  # Grow forest
  forest <- ranger.stab(data.bin.y, num.trees.step, alpha, mtry, num.trees, num.threads, replace,
                        sample.fraction, verbose, seed)
  paths <- forest$paths
  proba <- forest$proba
  num.trees <- forest$num.trees

  # path selection with p0
  if (!is.null(p0)) {
    selector <- proba > p0
    paths <- paths[selector]
    proba <- proba[selector]
    num.rule <- num.rule.max
  }

  # path post-treatment
  paths.select <- paths.filter(paths, proba, num.rule)
  paths <- paths.select$paths
  proba <- paths.select$proba

  # recover rule constraints from paths
  rules <- lapply(paths, get.rule, quantiles.emp = quantiles.emp, data.names = data.names)

  # rule outputs
  rules.out <- lapply(paths, get.rule.outputs, data.bin = data.bin, y = y)
  mean.out <- mean(y)

  # format paths
  paths <- lapply(paths, function(path){
    lapply(path, function(split){
      quantiles <- quantiles.emp[[split[1]]][[1]]
      quantiles.dedup <- unique(quantiles)
      is.cat <- if(length(quantiles.dedup)==2){all(quantiles.dedup == c(0,1))}else{FALSE}
      if (is.cat){
        split[2] <- split[3]
        split[3] <- '='
      }else{
        split.value <- quantiles.dedup[ceiling(split[2])]
        split[2] <- min(which(quantiles == split.value)) - 1
        split[3] <- if(split[3]==0){'L'}else{'R'}
      }
      return(split)
    })
  })

  return(list(rules = rules, rules.out = rules.out, proba = proba, paths = paths, num.trees = num.trees,
              mean = mean.out))

}

#'
#' Print the list of rules output by SIRUS.
#'
#' @title Print SIRUS
#' @param sirus.m A SIRUS model generated by sirus.fit.
#'
#' @return Formatted list of rules.
#' @export
#'
#' @examples
#' ## load sirus
#' require(sirus)
#'
#' ## prepare data
#' data <- iris
#' y <- rep(0, nrow(data))
#' y[data$Species == 'setosa'] = 1
#' data$Species <- NULL
#'
#' ## fit sirus
#' sirus.m <- sirus.fit(data, y)
#'
#' ## print sirus model
#' sirus.print(sirus.m)
#'
sirus.print <- function(sirus.m){

  # check sirus.m is a valid sirus model
  sirus.m.valid <- sirus.model.check(sirus.m)
  if (!sirus.m.valid){
    stop('Invalid sirus model.')
  }

  rules <- sirus.m$rules
  rules.out <- sirus.m$rules.out

  # format sirus output in a readable format
  if (length(rules) > 0){
    rules.print <- paste0(lapply(1:length(rules), function(ind) {
      rule <- rules[[ind]]
      rule.paste <- paste0(lapply(rule, function(split){
        split[3] <- signif(as.numeric(split[3]), digits = 3)
        paste0(split, collapse = ' ')
      }), collapse = ' & ')
      rule.out <- rules.out[[ind]]
      out.true <- signif(rule.out$outputs[1], digits = 3)
      out.false <- signif(rule.out$outputs[2], digits = 3)
      paste0(c('if', rule.paste, 'then', out.true, 'else', out.false), collapse = ' ')
    }))
    mean.print <- paste('Proportion of class 1:', signif(sirus.m$mean, digits = 3))
    rules.print <- c(mean.print, rules.print)
  }else{
    out <- signif(sirus.m$mean, digits = 3)
    rules.print <- paste0('Empty rule set. Constant output = ', out)
  }

  return(rules.print)
}

#'
#' Predictions of a SIRUS model for new observations.
#'
#' @title Predict
#' @param sirus.m A SIRUS model generated by sirus.fit.
#' @param data.test Testing data (dataframe of new observations).
#'
#' @return Predictions. A vector of the predicted probability of each new observation to be of class 1.
#' @export
#'
#' @examples
#' ## load sirus
#' require(sirus)
#'
#' ## prepare data
#' data <- iris
#' y <- rep(0, nrow(data))
#' y[data$Species == 'setosa'] = 1
#' data$Species <- NULL
#'
#' #' ## fit sirus
#' sirus.m <- sirus.fit(data, y)
#'
#' ## predict
#' predictions <- sirus.predict(sirus.m, data)
#'
sirus.predict <- function(sirus.m, data.test){

  # check sirus.m is a valid sirus model
  sirus.m.valid <- sirus.model.check(sirus.m)
  if (!sirus.m.valid){
    stop('Invalid sirus model.')
  }
  # check data.test
  data.check(data.test)

  rules <- sirus.m$rules
  rules.out <- sirus.m$rules.out

  if (length(rules) > 0){
    data.rule <- get.data.rule(data.test, rules, rules.out)
    pred <- apply(data.rule, 1, mean)
  }else{
    pred <- rep(sirus.m$mean, nrow(data.test))
  }

  return(pred)

}

#'
#' Estimation by cross-validation of the hyperparameter p0 used to select rules in sirus.fit.
#' For a robust estimation, it is recommanded to run multiple cross-validations, typically 30.
#'
#' @title Estimation of p0.
#' @param data Input dataframe, each row is an observation vector.
#' @param y Binary response variable taking 0 and 1 values.
#' @param nfold Number of folds in the cross-validation. Default is 10.
#' @param ncv Number of repetitions of the cross-validation. Default is 30.
#' @param num.rule.max Maximum number of rules of SIRUS model in the cross-validation grid. Default is 25.
#' @param q Number of quantiles used for node splitting in the forest construction. Default is 10.
#' @param num.trees.step Number of trees grown between two evaluations of the stopping criterion. Ignored if num.trees is provided.
#' @param alpha Parameter of the stopping criterion for the number of trees: stability has to reach 1 - alpha to stop the growing of the forest. Ignored if num.trees is provided.
#' @param mtry Number of variables to possibly split at each node. Default is the number of variables divided by 3.
#' @param verbose Boolean. If true, information messages are printed.
#' @param num.trees Number of trees grown in the forest. If NULL (recommanded), the number of trees is automatically set using a stability stopping criterion.
#' @param num.threads Number of threads used to grow the forest. Default is number of CPUs available.
#' @param replace Boolean. If true (default), sample with replacement.
#' @param sample.fraction Fraction of observations to sample. Default is 1 for sampling with replacement and 0.632 for sampling without replacement.
#' @param seed Random seed. Default is NULL, which generates the seed from R. Set to 0 to ignore the R seed.
#'
#' @return Optimal value of p0 with the elements
#'   \item{\code{p0}}{Optimal p0 value.}
#'   \item{\code{error}}{1 - AUC of the SIRUS model with the optimal p0.}
#'   \item{\code{num.rules}}{Number of rules of the SIRUS model with the optimal p0.}
#'   \item{\code{error.grid.p0}}{Table with the full cross-validation results for a fine grid of p0: number of rules, stability, error (1-AUC).}
#' @export
#'
#' @examples
#' ## load sirus
#' require(sirus)
#'
#' ## prepare data
#' data <- iris
#' y <- rep(0, nrow(data))
#' y[data$Species == 'setosa'] = 1
#' data$Species <- NULL
#'
#' ## run cv
#' cv.grid <- sirus.cv(data, y, nfold = 3, ncv = 2, num.trees = 100)
#'
sirus.cv <- function(data, y, nfold = 10, ncv = 30, num.rule.max = 25, q = 10,
                     num.trees.step = 1000, alpha = 0.05, mtry = NULL, num.trees = NULL,
                     num.threads = NULL, replace = TRUE, sample.fraction = NULL,
                     verbose = TRUE, seed = NULL){

  # check arguments
  data.y.check(data, y)
  #check nfold
  nfold.valid <- is.numeric(nfold)
  if (nfold.valid){nfold.valid <- (round(nfold) == nfold) & nfold >= 2 & nfold <= nrow(data)}
  if (!nfold.valid){
    stop('Invalid nfold. Number of cross-validation folds has to be an integer between 2 and the sample size.')
  }
  # check ncv
  ncv.valid <- is.numeric(ncv)
  if (ncv.valid){ncv.valid <- (round(ncv) == ncv) & ncv >= 1}
  if (!ncv.valid){
    stop('Invalid ncv. Number of cross-validations has to be an integer greater than 1.')
  }else{
    if (ncv == 1){
      warning('Warning ncv: It is recommanded to run multiple cross-validations for a robust estimation of p0.')
    }
  }
  # check sirus parameters
  sirus.param.check(data, num.rule.max, q, num.trees.step, alpha, mtry)
  # set default mtry
  if (is.null(mtry)){
    mtry.ratio <- 1/3
    mtry <- floor((ncol(data))*mtry.ratio)
  }else{
    mtry.ratio <- mtry/ncol(data)
  }
  # set sample fraction
  if (is.null(sample.fraction)){
    sample.fraction <- ifelse(replace, 1, 0.632)
  }
  ## Seed
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # p0.grid
  p0.grid <- exp(seq(log(mtry.ratio), log(1/1000), length.out = 500))

  # run cross-validation
  error.grids <- lapply(1:ncv, function(iter, p0.grid) {

    if (verbose == TRUE){
      print(paste0('Running cross-validation ', iter, '/', ncv, ' ...'))
    }

    # create k-folds
    ndata <- nrow(data)
    ind <- cut(seq(1,ndata), breaks = nfold, labels = F)
    ind <- sample(ind, size = ndata, replace = F)
    folds.ind <- lapply(1:nfold, function(fold, ind){
      test <- which(ind == fold)
      train <- setdiff(1:ndata, test)
      return(list(train = train, test = test))
    }, ind  = ind)

    # fit sirus for each fold and compute prediction for all rule selections
    pred.cv <- lapply(1:nfold, function(fold){

      data.train <- data[folds.ind[[fold]]$train,]
      y.train <- y[folds.ind[[fold]]$train]
      sirus.cv <- sirus.fit(data.train, y.train, num.rule = (num.rule.max + 10), p0 = NULL, q = q,
                            num.trees.step = num.trees.step, alpha = alpha, mtry = mtry, num.trees = num.trees,
                            num.threads = num.threads, replace = replace, sample.fraction = sample.fraction,
                            verbose = FALSE, seed = seed)
      data.test <- data[folds.ind[[fold]]$test,]
      data.rule.cv <- get.data.rule(data.test, sirus.cv$rules, sirus.cv$rules.out)
      pred.df <- sapply(1:ncol(data.rule.cv), function(ind){
        apply(data.rule.cv[,1:ind,drop=F], 1, mean)
      })
      pred.df <- cbind(rep(sirus.cv$mean, nrow(pred.df)), pred.df)
      proba.cv <- sirus.cv$proba
      paths.cv <- sirus.cv$paths

      return(list(pred.df = pred.df, proba = proba.cv, paths = paths.cv))

    })

    # compute prediction error (1-AUC)
    proba.cv <- lapply(pred.cv, function(pred){pred$proba})
    pred.df <- lapply(pred.cv, function(pred){
      pred$pred.df
    })
    y.test <- unlist(lapply(1:nfold, function(fold){
      y[folds.ind[[fold]]$test]
    }))
    folds.ncol <- sapply(1:nfold, function(fold){
      proba <- proba.cv[[fold]]
      proba <- c(proba, 0)
      fold.ncol <- rep(0, length(p0.grid))
      for (ind in 2:(min(length(proba), (num.rule.max + 10) + 1))){
        fold.ncol[proba[ind] <= p0.grid & proba[ind - 1] > p0.grid] <- ind - 1
      }
      fold.ncol
    })
    # TODO: dedup identical computations to speed-up
    error.grid <- unlist(lapply(1:length(p0.grid), function(ind){
      pred <- unlist(lapply(1:nfold, function(fold){
        pred.df[[fold]][,folds.ncol[ind,fold]+1]
      }))
      1 - as.numeric(performance(prediction(pred, y.test), "auc")@y.values)
    }))

    # compute stability metric
    stab.df <- lapply(1:(nfold - 1), function(fold1){
      pred.cv1 <- pred.cv[[fold1]]
      sapply((fold1 + 1):nfold, function(fold2){
        pred.cv2 <- pred.cv[[fold2]]
        sapply(p0.grid, function(p0){
          paths1 <- pred.cv1$paths[pred.cv1$proba > p0]
          paths2 <- pred.cv2$paths[pred.cv2$proba > p0]
          len <- length(paths1) + length(paths2)
          if (len > 0){
            2*length(intersect(paths1, paths2))/len
          }else{
            1
          }
        })
      })
    })
    stab.df <- do.call('cbind', stab.df)
    stab.grid <- apply(stab.df, 1, mean)

    # compute number of rules
    num.rules <- apply(folds.ncol, 1, mean)

    return(list(error.grid = error.grid, stab.grid = stab.grid, num.rules = num.rules))

  }, p0.grid = p0.grid)

  # aggregate results
  error.df <- sapply(error.grids, function(x){x$error.grid})
  error.mean <- apply(error.df, 1, mean)
  if (ncv > 1){
    error.sd <- apply(error.df, 1, sd)/sqrt(ncv)
  }else{
    error.sd <- rep(0, length(p0.grid))
  }
  stab.df <- sapply(error.grids, function(x){x$stab.grid})
  stab.mean <- apply(stab.df, 1, mean)
  if (ncv > 1){
    stab.sd <- apply(stab.df, 1, sd)/sqrt(ncv)
  }else{
    stab.sd <- rep(0, length(p0.grid))
  }
  num.rules.df <- sapply(error.grids, function(x){x$num.rules})
  num.rules.mean <- apply(num.rules.df, 1, mean)
  if (ncv > 1){
    num.rules.sd <- apply(num.rules.df, 1, sd)/sqrt(ncv)
  }else{
    num.rules.sd <- rep(0, length(p0.grid))
  }
  error.grid.p0 <- as.data.frame(cbind(p0.grid, num.rules.mean, stab.mean, error.mean,
                                       num.rules.sd, stab.sd, error.sd))
  ind.1 <- which.min(abs(1 - num.rules.mean))[1]
  ind.max <- which.min(abs(min(num.rule.max, max(num.rules.mean)) - num.rules.mean))
  error.grid.p0 <- error.grid.p0[ind.1:ind.max,]

  ind.min <- which.min(error.mean[1:ind.max])
  ind.opt <- min(which(error.mean <= error.mean[ind.min] + 2*error.sd[ind.min]))
  p0.opt <- p0.grid[ind.opt]
  num.rules.opt <- num.rules.mean[ind.opt]
  auc.opt <- error.mean[ind.opt]
  stab.opt <- stab.mean[ind.opt]

  return(list(p0 = p0.opt, error = auc.opt, stab = stab.opt, num.rules = num.rules.opt,
              error.grid.p0 = error.grid.p0))

}

#'
#' Plot SIRUS cross-validation path: 1-AUC and stability versus the number of rules when p0 varies.
#'
#' @title Plot SIRUS cross-validation.
#' @param sirus.cv.grid Cross-validation results returned by sirus.cv.
#' @param num.rule.max Upper limit on the number of rules for the x-axis. Default is 25.
#'
#' @return Plots of cross-validation results.
#'   \item{\code{error}}{plot of 1-AUC vs number of rules (ggplot object).}
#'   \item{\code{stab}}{plot of stability vs number of rules (ggplot object).}
#' @export
#'
#' @examples
#' ## load sirus
#' require(sirus)
#'
#' ## prepare data
#' data <- iris
#' y <- rep(0, nrow(data))
#' y[data$Species == 'setosa'] = 1
#' data$Species <- NULL
#'
#' ## run cv
#' cv.grid <- sirus.cv(data, y, nfold = 3, ncv = 2, num.trees = 100)
#'
#' ## plot cv result
#' plot.error <- sirus.plot.cv(cv.grid)$error
#' plot(plot.error)
#'
sirus.plot.cv <- function(sirus.cv.grid, num.rule.max = 25){

  # filter cv grid
  error.grid.p0 <- sirus.cv.grid$error.grid.p0
  num.rule.max <- min(num.rule.max, max(error.grid.p0$num.rules.mean))
  grid.index <- sapply(1:num.rule.max, function(ind){which.min(abs(ind - sirus.cv.grid$error.grid.p0$num.rules.mean))[1]})
  error.grid.p0 <- sirus.cv.grid$error.grid.p0[grid.index,]

  # declare variables
  num.rules.mean <- NULL
  error.mean <- NULL
  error.sd <- NULL
  stab.mean <- NULL
  stab.sd <- NULL

  # plot 1 - AUC vs number of rules
  tag.y <- (min(error.grid.p0$error.mean) + max(error.grid.p0$error.mean))/2
  plot.error <- ggplot(error.grid.p0, aes(x = num.rules.mean, y = error.mean)) +
    geom_line(size = 0.8) + geom_point(size=1) +
    geom_errorbar(aes(ymin=error.mean - 2*error.sd,
                      ymax=error.mean + 2*error.sd),  width=0.7, size=0.5)+
    geom_hline(yintercept = sirus.cv.grid$error,
               linetype = 'dashed', color = 'blue', size = 0.7) +
    geom_vline(xintercept = sirus.cv.grid$num.rules,
               linetype = 'dashed', color = 'blue', size = 0.7) +
    geom_text(aes(sirus.cv.grid$num.rules, tag.y, label = 'Optimal p0'), angle = '90',
              vjust=1, color='blue', size = 7) +
    geom_text(aes(num.rule.max - 5, sirus.cv.grid$error, label = paste0('SIRUS - ', round(sirus.cv.grid$error,2))),
              vjust=-1, color='blue', size = 7) +
    xlab('Number of rules') +
    ylab('1-AUC') +
    theme_classic() +
    theme(axis.line.x = element_line(colour = 'black'),
          axis.line.y = element_line(colour = 'black'),
          text = element_text(size=20),
          plot.title = element_text(hjust = 0.5, size=18, face="italic"))

  # plot stability vs number of rules
  plot.stab <- ggplot(error.grid.p0, aes(x = num.rules.mean, y = stab.mean)) +
    geom_line(size = 0.8) + geom_point() +
    geom_errorbar(aes(ymin=stab.mean - stab.sd,
                      ymax=stab.mean + stab.sd),  width=0.7, size=0.5)+
    geom_vline(xintercept = sirus.cv.grid$num.rules,
               linetype = 'dashed', color = 'blue', size = 0.7) +
    geom_text(aes(sirus.cv.grid$num.rules, mean(error.grid.p0$stab.mean), label = 'Optimal p0'), angle = '90',
              vjust=1, color='blue', size = 7) +
    xlab('Number of rules') +
    ylab('Stability') +
    theme_classic() + ylim(c(0.9*min(error.grid.p0$stab.mean), 1.1*max(error.grid.p0$stab.mean)))
    theme(axis.line.x = element_line(colour = 'black'),
        axis.line.y = element_line(colour = 'black'),
        text = element_text(size=20),
        plot.title = element_text(hjust = 0.5, size=18, face="italic"))

  return(list(error = plot.error, stability = plot.stab))

}
