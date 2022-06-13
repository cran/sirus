#'
#' Fit SIRUS for a given number of rules (10 by default) or a given \code{p0}. \cr
#' SIRUS is a regression and classification algorithm, based on random forests (Breiman, 2001), that takes the form of a short list of rules.
#' SIRUS combines the simplicity of rule algorithms or decision trees with an accuracy close to random forests.
#' More importantly, the rule selection is stable with respect to data perturbation.
#' SIRUS for classification is defined in (Benard et al. 2021a), and the extension to regression is provided in (Benard et al. 2021b).
#' 
#' If the output \code{y} takes only 0 and 1 values, a classification model is fit, otherwise a regression model is fit.
#' SIRUS algorithm proceeds the following steps:
#' \enumerate{
#'   \item Discretize data
#'   \item Fit a random forest
#'   \item Extract rules from tree nodes
#'   \item Select the most frequent rules (which occur in at least a fraction p0 of the trees)
#'   \item Filter rules to remove linear dependence between them
#'   \item Aggregate the selected rules
#'      \itemize{
#'         \item Classification: rules are averaged
#'         \item Regression: rules are linearly combined via a ridge regression (constrained to have all coefficients positive)
#'      }
#' }
#' The hyperparameter \code{p0} can be tuned using \code{\link{sirus.cv}} to set the optimal number of rules. \cr
#' The number of trees is automatically set with a stopping criterion based on stability: 
#' the forest growing is stopped when the number of trees is high enough to ensure that 95\% of the rules in average are identical over two runs of SIRUS on the provided dataset. \cr
#' Data is discretized depending on variable types: numerical variables are binned using \code{q}-quantiles, 
#' categorical variables are transformed in ordered variables as in \code{\link[ranger]{ranger}} (standard method to handle categorical variables in trees),
#' while discrete variables (numerical variables with less than \code{discrete.limit} distinct values) are left untouched.
#' Notice that categorical variables with a high number of categories should be discarded or transformed, as SIRUS is likely to identify associated irrelevant rules.
#'
#' @title Fit SIRUS.
#'
#' @param data Input dataframe, each row is an observation vector. Each column is an input variable and is numeric or factor.
#' @param y Numeric response variable. For classification, \code{y} takes only 0 and 1 values.
#' @param type 'reg' for regression, 'classif' for classification and 'auto' for automatic detection (classification if \code{y} takes only 0 and 1 values).
#' @param num.rule Number of rules in SIRUS model. Default is 10. Ignored if a \code{p0} value is provided. For regression, the effective number of rules can be smaller than \code{num.rule} because of null coefficients in the final linear aggregation of the rules.
#' @param p0 Selection threshold on the frequency of appearance of a path in the forest to set the number of rules. Default is NULL and \code{num.rule} is used to select rules. \code{\link{sirus.cv}} provides the optimal \code{p0} by cross-validation.  
#' @param num.rule.max Maximum number of rules in SIRUS model. Ignored if \code{num.rule} is provided.
#' @param q Number of quantiles used for node splitting in the forest construction. Default and recommended value is 10.
#' @param discrete.limit Maximum number of distinct values for a variable to be considered discrete. If higher, variable is continuous.
#' @param num.trees.step Number of trees grown between two evaluations of the stopping criterion. Ignored if \code{num.trees} is provided.
#' @param alpha Parameter of the stopping criterion for the number of trees: stability has to reach 1-\code{alpha} to stop the growing of the forest. Ignored if \code{num.trees} is provided. Default value is 0.05.
#' @param mtry Number of variables to possibly split at each node. Default is the number of variables divided by 3.
#' @param max.depth Maximal tree depth. Default and recommended value is 2.
#' @param num.trees Number of trees grown in the forest. Default is NULL. If NULL (recommended), the number of trees is automatically set using a stability based stopping criterion.
#' @param num.threads Number of threads used to grow the forest. Default is number of CPUs available.
#' @param replace Boolean. If true (default), sample with replacement.
#' @param sample.fraction Fraction of observations to sample. Default is 1 for sampling with replacement and 0.632 for sampling without replacement.
#' @param verbose Boolean. If true, information messages are printed.
#' @param seed Random seed. Default is NULL, which generates the seed from R. Set to 0 to ignore the R seed.
#' 
#' @return SIRUS model with elements
#'   \item{\code{rules}}{List of rules in SIRUS model.}
#'   \item{\code{rules.out}}{List of rule outputs. \code{rule.out}: the output mean whether the rule is satisfied or not. \code{supp.size}: the number of points inside and outside the rule.}
#'   \item{\code{proba}}{Frequency of occurence of paths in the forest.}
#'   \item{\code{paths}}{List of selected paths (symbolic representation with quantile order for continuous variables).}
#'   \item{\code{rule.weights}}{Vector of positive or null coefficients assigned to each rule for the linear aggregation (1/number of rules for classification).}
#'   \item{\code{rule.glm}}{Fitted glmnet object for regression (linear rule aggregation with ridge penalty).}
#'   \item{\code{type}}{Type of SIRUS model: 'reg' for regression, 'classif' for classification.} 
#'   \item{\code{num.trees}}{Number of trees used to build SIRUS.}
#'   \item{\code{data.names}}{Names of input variables.}
#'   \item{\code{mean}}{Mean output over the full training data. Default model output if no rule is selected.}
#'   \item{\code{bins}}{List of type and possible split values for all input variables.}
#' @export
#'
#' @examples
#' ## load SIRUS
#' require(sirus)
#'
#' ## prepare data
#' data <- iris
#' y <- rep(0, nrow(data))
#' y[data$Species == 'setosa'] = 1
#' data$Species <- NULL
#'
#' ## fit SIRUS
#' sirus.m <- sirus.fit(data, y)
#'
#' @references
#' \itemize{
#'   \item Benard, C., Biau, G., Da Veiga, S. & Scornet, E. (2021a). SIRUS: Stable and Interpretable RUle Set for Classification. Electronic Journal of Statistics, 15:427-505. \doi{10.1214/20-EJS1792}.
#'   \item Benard, C., Biau, G., Da Veiga, S. & Scornet, E. (2021b). Interpretable Random Forests via Rule Extraction. Proceedings of The 24th International Conference on Artificial Intelligence and Statistics, PMLR 130:937-945. \url{http://proceedings.mlr.press/v130/benard21a}. 
#'   \item Breiman, L. (2001). Random forests. Machine learning, 45, 5-32.
#'   \item Wright, M. N. & Ziegler, A. (2017). ranger: A fast implementation of random forests for high dimensional data in C++ and R. J Stat Softw 77:1-17. \doi{10.18637/jss.v077.i01}.
#' }
#' @encoding UTF-8
#' @useDynLib sirus
#' @importFrom Rcpp evalCpp
#' @import stats
#' @import utils
#' @importFrom Matrix Matrix
#' @importFrom Matrix rankMatrix
#' @import ROCR
#' @import ggplot2
#' @import glmnet
sirus.fit <- function(data, y, type = 'auto', num.rule = 10, p0 = NULL, num.rule.max = 25, q = 10, 
                      discrete.limit = 10, num.trees.step = 1000, alpha = 0.05, mtry = NULL,
                      max.depth = 2, num.trees = NULL, num.threads = NULL, replace = TRUE,
                      sample.fraction = ifelse(replace, 1, 0.632), verbose = TRUE, seed = NULL) {

  # Check arguments
  # check data type
  data.y.check(data, y)
  # check SIRUS parameters
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
    mtry <- max(1, floor((ncol(data))*mtry.ratio))
  }
  # set type
  type.valid <- type %in% c('auto', 'reg', 'classif')
  if (type.valid){
    if (type == 'reg'){
      if (length(unique(y)) == 2){
        warning('Warning type: y takes only two distinct values, classification is more appropriate.')
      }
    }
    if (type == 'classif'){
      if (!all(unique(y) %in% c(0, 1))){
        stop('Invalid type: for classification, y has to take only 0 and 1 values.')
      }
    }
    if (type == 'auto'){
      if (all(unique(y) %in% c(0, 1))){type = 'classif'}else{type = 'reg'}
    }
  }else{
    stop('Invalid type. Type should be auto, reg (for regression) or classif (for classification).')
  }

  # Data info
  data.names <- colnames(data)
  mean.out <- mean(y)
  
  # Data binning
  bins.list <- lapply(data, get.bins, y = y, q = q, discrete.limit = discrete.limit)
  num.cat.valid <- sapply(bins.list, function(bins){
    if (bins$type =='categorical'){length(bins$levels) < 0.05*nrow(data)}else{TRUE}
  })
  if (!all(num.cat.valid)){
    names.warning <- paste0(data.names[!num.cat.valid], collapse = ', ')
    warning(paste0('Some categorical variables have a high number of categories. 
    SIRUS is likely to identify irrelevant rules, overfit, and to have long running times.
    It is recommended to discard or transform these variables: ', names.warning, '.'))
  }
  data.bin <- as.data.frame(sapply(1:ncol(data), function(j){
                              binarize.X(X = data[,j], bins = bins.list[[j]], q = q)}))
  data.bin.y <- cbind(data.bin, y)

  if (nrow(data) > 5){
    # Grow forest
    forest <- ranger.stab(data.bin.y, num.trees.step, alpha, mtry, max.depth, num.trees,
                          num.threads, replace, sample.fraction, verbose, seed)
    paths <- forest$paths[-1]
    proba <- forest$proba[-1]
    num.trees <- forest$num.trees
  }else{
    warning('Minimum sample size to run SIRUS is 5.')
    paths <- list()
    proba <- 1
    num.trees <- 0
  }

  if (length(paths) > 0){
    
    # path selection with p0
    if (!is.null(p0)) {
      selector <- proba > p0
      paths <- paths[selector]
      proba <- proba[selector]
      num.rule <- num.rule.max
    }
    
    # path post-treatment
    paths <- lapply(paths, function(path){
      lapply(path, function(split){
        if (split[2] == round(split[2])){
          split[2] <- split[2] - 0.5
        }
        split
      })
    })
    if (max.depth <= 2){
      paths.ftr <- paths.filter.2(paths, proba, num.rule)
    }else{
      paths.ftr <- paths.filter.d(paths, proba, num.rule, data.bin)
    }
    paths <- paths.ftr$paths
    proba <- paths.ftr$proba
    
    # format paths
    paths <- lapply(paths, format.path, bins.list = bins.list)
  
    # build rules from paths
    rules <- lapply(paths, get.rule, bins.list = bins.list, data.names = data.names)
    data.rule.supp <- get.rule.support(data, rules)
    rules.out <- get.rule.outputs(data.rule.supp, y)
    
    # symbolic paths
    paths <- lapply(paths, function(path){lapply(path, function(split){
                if (bins.list[[as.numeric(split[1])]]$type %in% c('continuous', 'discrete')){
                  split <- split[1:3]}; return(split)})
              })
  
    # for regreesion: fit rule coefficients
    if (type == 'reg' & length(paths) > 1){
      data.rule <- sapply(1:length(paths), function(j){
        X <- data.rule.supp[, j]
        X[X == 1] <- rules.out[[j]]$outputs[1]
        X[X == 0] <- rules.out[[j]]$outputs[2]
        X
      })
      lambda <- cv.glmnet(data.rule, y, lower.limits = rep(0, ncol(data.rule)),
                          alpha = 0, standardize = F)$lambda.min
      rule.glm <- glmnet(data.rule, y, lower.limits = rep(0, ncol(data.rule)),
                          lambda = lambda, alpha = 0, standardize = F)
      rule.weights <- as.vector(rule.glm$beta)
    }else{
      rule.glm <- NULL
      rule.weights <- rep(1/length(rules), length(rules))
    }
    
  }else{
    
    rules <- list()
    rules.out <- list()
    rule.glm <- NULL
    rule.weights <- 1
    
  }

  return(list(rules = rules, rules.out = rules.out, proba = proba, paths = paths,
              rule.weights = rule.weights, rule.glm = rule.glm, type = type, 
              num.trees = num.trees, data.names = data.names, mean = mean.out, bins = bins.list))

}

#'
#' Print the list of rules output by SIRUS.
#'
#' @title Print SIRUS.
#'
#' @param sirus.m A SIRUS model generated by \code{\link{sirus.fit}}.
#' @param digits Number of significant digits for numerical values. Default value is 3.
#'
#' @return Formatted list of rules.
#' @export
#'
#' @examples
#' ## load SIRUS
#' require(sirus)
#'
#' ## prepare data
#' data <- iris
#' y <- rep(0, nrow(data))
#' y[data$Species == 'setosa'] = 1
#' data$Species <- NULL
#'
#' ## fit SIRUS
#' sirus.m <- sirus.fit(data, y)
#'
#' ## print sirus model
#' sirus.print(sirus.m)
#'
sirus.print <- function(sirus.m, digits = 3){

  # check sirus.m is a valid SIRUS model
  sirus.m.valid <- sirus.model.check(sirus.m)
  if (!sirus.m.valid){
    stop('Invalid SIRUS model.')
  }
  digits.valid <- is.numeric(digits) & digits > 0
  if (!digits.valid){
    stop('Invalid number of digits.')
  }

  rules <- sirus.m$rules
  rules.out <- sirus.m$rules.out
  rule.weights <- sapply(sirus.m$rule.weights, signif, digits = digits)

  # format SIRUS output in a readable format
  if (length(rules) > 0){
    rules.print <- paste0(lapply(1:length(rules), function(ind) {
      rule <- rules[[ind]]
      rule.paste <- paste0(lapply(rule, function(split){
        if (split[2] %in% c('<', '>=')){
          split[3] <- signif(as.numeric(split[3]), digits = digits)
          split.paste <- paste0(split, collapse = ' ')
        }
        if (split[2] == '='){
          split.paste <- paste0(split[1], ' in {', paste0(split[3:length(split)], collapse = ', '), '}')
        }
        split.paste
      }), collapse = ' & ')
      rule.out <- rules.out[[ind]]
      out.true <- signif(rule.out$outputs[1], digits = digits)
      out.false <- signif(rule.out$outputs[2], digits = digits)
      size.true <- rule.out$supp.size[1]
      size.false <- rule.out$supp.size[2]
      weight <- signif(rule.weights[ind], digits = digits)
      paste0(c('if ', rule.paste, ' then ', out.true, ' (n=', size.true, ') else ', out.false, ' (n=', size.false, ')'), collapse = '')
    }))
    
    if (sirus.m$type == 'classif'){
      mean.print <- paste('Proportion of class 1 =', signif(sirus.m$mean, digits = digits),
                          '- Sample size n =', sum(rules.out[[1]]$supp.size))
      rules.print <- c(mean.print, rules.print)
    }else{
      rules.print <- rules.print[rule.weights > 0]
      rule.weights <- rule.weights[rule.weights > 0]
      mean.print <- paste0('Mean of output y = ', signif(sirus.m$mean, digits = digits),
                          ' - Sample size n = ', sum(rules.out[[1]]$supp.size))
      rules.print <- c(mean.print, rules.print)
      if (length(rules) > 1){intercept <- sirus.m$rule.glm$a0}else{intercept <- 0}
      intercept <- paste('Intercept =', signif(intercept, digits = digits))
      rule.weights <- c(intercept, rule.weights)
      rules.print <- cbind(rule.weights, rules.print)
      colnames(rules.print) <- c('Weights', 'Rules')
    }
    
  }else{
    out <- signif(sirus.m$mean, digits = digits)
    rules.print <- paste0('Empty rule set. Constant output = ', out)
  }

  return(rules.print)
}

#'
#' Compute SIRUS predictions for new observations.
#'
#' @title Predict.
#' @param sirus.m A SIRUS model generated by \code{\link{sirus.fit}}.
#' @param data.test Testing data (dataframe of new observations).
#'
#' @return Predictions. For classification, vector of the predicted probability of each new observation to be of class 1.
#' @export
#'
#' @examples
#' ## load SIRUS
#' require(sirus)
#'
#' ## prepare data
#' data <- iris
#' y <- rep(0, nrow(data))
#' y[data$Species == 'setosa'] = 1
#' data$Species <- NULL
#'
#' #' ## fit SIRUS
#' sirus.m <- sirus.fit(data, y)
#'
#' ## predict
#' predictions <- sirus.predict(sirus.m, data)
#'
sirus.predict <- function(sirus.m, data.test){

  # check sirus.m is a valid SIRUS model
  sirus.m.valid <- sirus.model.check(sirus.m)
  if (!sirus.m.valid){
    stop('Invalid SIRUS model.')
  }
  # check data.test
  data.check(data.test)
  # check names
  if (!is.null(colnames(data.test))){
    names.valid <- all(colnames(data.test) == sirus.m$data.names)
  }else{names.valid <- FALSE}
  if (!names.valid){
    stop('Invalid variable names for data.test.')
  }

  rules <- sirus.m$rules
  rules.out <- sirus.m$rules.out

  if (length(rules) > 0){
    data.rule <- get.data.rule(data.test, rules, rules.out)
    if (sirus.m$type == 'classif' | length(rules) == 1){
      pred <- apply(data.rule, 1, mean)
    }else{
      pred <- as.vector(predict(sirus.m$rule.glm, data.rule))
    }
  }else{
    pred <- rep(sirus.m$mean, nrow(data.test))
  }

  return(pred)

}

#'
#' Estimate the optimal hyperparameter \code{p0} used to select rules in \code{\link{sirus.fit}} using cross-validation (Benard et al. 2021a, 2021b).
#' 
#' For a robust estimation of \code{p0}, it is recommended to run multiple cross-validations (typically \code{ncv} = 10).
#' Two optimal values of \code{p0} are provided: \code{p0.pred} (Benard et al. 2021a) and \code{p0.stab} (Benard et al. 2021b), defined such that \code{p0.pred} minimizes the error, and \code{p0.stab} finds a tradeoff between error and stability.
#' Error is 1-AUC for classification and the unexplained variance for regression.
#' Stability is the average proportion of rules shared by two SIRUS models fit on two distinct folds of the cross-validation.
#'
#' @title Estimate p0.
#'
#' @param data Input dataframe, each row is an observation vector. Each column is an input variable and is numeric or factor.
#' @param y Numeric response variable. For classification, \code{y} takes only 0 and 1 values.
#' @param type 'reg' for regression, 'classif' for classification and 'auto' for automatic detection (classification if \code{y} takes only 0 and 1 values).
#' @param nfold Number of folds in the cross-validation. Default is 10.
#' @param ncv Number of repetitions of the cross-validation. Default is 10 for a robust estimation of \code{p0}.
#' @param num.rule.max Maximum number of rules of SIRUS model in the cross-validation grid. Default is 25.
#' @param q Number of quantiles used for node splitting in the forest construction. Default and recommended value is 10.
#' @param discrete.limit Maximum number of distinct values for a variable to be considered discrete. If higher, variable is continuous.
#' @param num.trees.step Number of trees grown between two evaluations of the stopping criterion. Ignored if \code{num.trees} is provided.
#' @param alpha Parameter of the stopping criterion for the number of trees: stability has to reach 1-\code{alpha} to stop the growing of the forest. Ignored if \code{num.trees} is provided. Default value is 0.05.
#' @param mtry Number of variables to possibly split at each node. Default is the number of variables divided by 3.
#' @param max.depth Maximal tree depth. Default and recommended value is 2.
#' @param num.trees Number of trees grown in the forest. If NULL (recommended), the number of trees is automatically set using a stability stopping criterion.
#' @param num.threads Number of threads used to grow the forest. Default is number of CPUs available.
#' @param replace Boolean. If true (default), sample with replacement.
#' @param sample.fraction Fraction of observations to sample. Default is 1 for sampling with replacement and 0.632 for sampling without replacement.
#' @param verbose Boolean. If true, information messages are printed.
#' @param seed Random seed. Default is NULL, which generates the seed from R. Set to 0 to ignore the R seed.
#'
#' @return Optimal value of \code{p0} with the elements
#'   \item{\code{p0.pred}}{Optimal \code{p0} value to minimize model error (recommended for classification).}
#'   \item{\code{p0.stab}}{Optimal \code{p0} value for a tradeoff between error and stability (recommended for regression).}
#'   \item{\code{error.grid.p0}}{Table with the full cross-validation results for a fine grid of \code{p0}: number of rules, stability, and error. 
#'                               The last three columns of the table are the standard deviations of the metrics across the \code{ncv} repetitions of the cross-validation.
#'                               See details for the definitions of the error and stability metrics.}
#'   \item{\code{type}}{'reg' for regression, 'classif' for classification.}
#' @export
#'
#' @examples
#' ## load SIRUS
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
#' @references
#' \itemize{
#'   \item Benard, C., Biau, G., Da Veiga, S. & Scornet, E. (2021a). SIRUS: Stable and Interpretable RUle Set for Classification. Electronic Journal of Statistics, 15:427-505. \doi{10.1214/20-EJS1792}.
#'   \item Benard, C., Biau, G., Da Veiga, S. & Scornet, E. (2021b). Interpretable Random Forests via Rule Extraction. Proceedings of The 24th International Conference on Artificial Intelligence and Statistics, PMLR 130:937-945. \url{http://proceedings.mlr.press/v130/benard21a}. 
#' }
sirus.cv <- function(data, y, type = 'auto', nfold = 10, ncv = 10, num.rule.max = 25, q = 10,
                     discrete.limit = 10, num.trees.step = 1000, alpha = 0.05, mtry = NULL,
                     max.depth = 2, num.trees = NULL, num.threads = NULL, replace = TRUE, 
                     sample.fraction = NULL, verbose = TRUE, seed = NULL){

  # check arguments
  data.y.check(data, y)
  #check nfold
  nfold.valid <- is.numeric(nfold)
  if (nfold.valid){nfold.valid <- (round(nfold) == nfold) & nfold >= 2}
  if (!nfold.valid){
    stop('Invalid nfold. Number of cross-validation folds has to be an integer greater than 2.')
  }else{
    if (nfold > nrow(data)){
      nfold <- nrow(data)
      warning(paste0('Warning nfold: nfold is greater than the sample size (=', nrow(data), '),
                      and is then automatically set to ', nfold, '.'))
    }
  }
  # set type
  type.valid <- type %in% c('auto', 'reg', 'classif')
  if (type.valid){
    if (type == 'reg'){
      if (length(unique(y)) == 2){
        warning('Warning type: y takes only two distinct values, classification is more appropriate.')
      }
    }
    if (type == 'classif'){
      if (!all(unique(y) %in% c(0, 1))){
        stop('Invalid type: for classification, y has to take only 0 and 1 values.')
      }
    }
    if (type == 'auto'){
      if (all(unique(y) %in% c(0, 1))){type = 'classif'}else{type = 'reg'}
    }
  }else{
    stop('Invalid type. Type should be auto, reg (for regression) or classif (for classification).')
  }
  # check ncv
  ncv.valid <- is.numeric(ncv)
  if (ncv.valid){ncv.valid <- (round(ncv) == ncv) & ncv >= 1}
  if (!ncv.valid){
    stop('Invalid ncv. Number of cross-validations has to be an integer greater than 1.')
  }else{
    if (ncv == 1){
      warning('Warning ncv: It is recommended to run multiple cross-validations for a robust estimation of p0.')
    }
  }
  # check SIRUS parameters
  sirus.param.check(data, num.rule.max, q, num.trees.step, alpha, mtry)
  # set default mtry
  if (is.null(mtry)){
    mtry.ratio <- 1/3
    mtry <- max(1, floor((ncol(data))*mtry.ratio))
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
    
    # max num rule
    if (type == 'classif'){
      num.rule.cv <- num.rule.max + 10
    }else{
      num.rule.cv <- 2*num.rule.max
    }

    # fit SIRUS for each fold and compute prediction for all rule selections
    pred.cv <- lapply(1:nfold, function(fold){

      data.train <- data[folds.ind[[fold]]$train, , drop = F]
      y.train <- y[folds.ind[[fold]]$train]
      sirus.cv <- sirus.fit(data.train, y.train, type = type, num.rule = num.rule.cv, p0 = NULL, q = q,
                            num.trees.step = num.trees.step, alpha = alpha, mtry = mtry, max.depth = max.depth,
                            num.trees = num.trees, num.threads = num.threads, replace = replace, 
                            sample.fraction = sample.fraction, verbose = FALSE, seed = seed)
      data.test <- data[folds.ind[[fold]]$test, , drop = F]
      if (length(sirus.cv$rules) > 0){
        data.rule.test <- get.data.rule(data.test, sirus.cv$rules, sirus.cv$rules.out)
        if (type == 'reg'){
          data.rule.train <- get.data.rule(data.train, sirus.cv$rules, sirus.cv$rules.out)
        }
      }else{
        data.rule.test <- as.data.frame(rep(sirus.cv$mean, nrow(data.test)))
        if (type == 'reg'){
          data.rule.train <- as.data.frame(rep(sirus.cv$mean, nrow(data.train)))
        }
      }
      pred.df <- lapply(1:ncol(data.rule.test), function(ind){
        if (type == 'classif'){
          pred <- apply(data.rule.test[,1:ind,drop=F], 1, mean)
          beta <- rep(1/ind, ind)
        }else{
          if (ind > 1){
            # Exception for last simu
            lambda <- cv.glmnet(data.rule.train[,1:ind,drop=F], y.train, lower.limits = rep(0, ind),
                                alpha = 0, standardize = F)$lambda.min
            sirus.glm <- glmnet(data.rule.train[,1:ind,drop=F], y.train, lower.limits = rep(0, ind),
                                lambda = lambda, alpha = 0, standardize = F)
            beta.tmp <- rep(0, num.rule.cv)
            beta.tmp[1:ind][as.vector(sirus.glm$beta) > 0] <- 1
            pred <- predict(sirus.glm, data.rule.test[,1:ind,drop=F])
            beta <- beta.tmp
          }else{
            beta.tmp <- c(1, rep(0, num.rule.cv - 1))
            pred <- data.rule.test[,1,drop=F]
            beta <- beta.tmp
          }
        }
        list(pred = pred, beta = beta)
      })
      beta.df <- lapply(pred.df, function(pred.ind){
        pred.ind$beta
      })
      pred.df <- matrix(sapply(pred.df, function(pred.ind){
        pred.ind$pred
      }), nrow = nrow(data.rule.test))
      pred.df <- cbind(rep(sirus.cv$mean, nrow(pred.df)), pred.df)
      proba.cv <- sirus.cv$proba
      paths.cv <- sirus.cv$paths

      return(list(pred.df = pred.df, beta.df = beta.df, proba = proba.cv, paths = paths.cv))

    })

    # compute prediction error (1-AUC/R2)
    proba.cv <- lapply(pred.cv, function(pred){pred$proba})
    pred.df <- lapply(pred.cv, function(pred){pred$pred.df})
    y.test <- unlist(lapply(1:nfold, function(fold){
      y[folds.ind[[fold]]$test]
    }))
    folds.ncol <- sapply(1:nfold, function(fold){
      proba <- proba.cv[[fold]]
      proba <- c(proba, 0)
      fold.ncol <- rep(0, length(p0.grid))
      for (ind in 2:(min(length(proba), num.rule.cv + 1))){
        fold.ncol[proba[ind] <= p0.grid & proba[ind - 1] > p0.grid] <- ind - 1
      }
      fold.ncol
    })
    # TODO: dedup identical computations to speed-up
    error.grid <- unlist(lapply(1:length(p0.grid), function(ind){
      pred <- unlist(lapply(1:nfold, function(fold){
        pred.df[[fold]][,folds.ncol[ind,fold]+1]
      }))
      if (type == 'classif'){
        if (length(unique(y.test)) == 1){
          0.5
        }else{
          1 - as.numeric(performance(prediction(pred, y.test), "auc")@y.values)
        }
      }else{
        sum((pred - y.test)^2)/sum((y.test - mean(y.test))^2)
      }
    }))

    # compute stability metric & number of rules
    stab.df <- lapply(1:(nfold - 1), function(fold1){
      pred.cv1 <- pred.cv[[fold1]]
      lapply((fold1 + 1):nfold, function(fold2){
        pred.cv2 <- pred.cv[[fold2]]
        sapply(p0.grid, function(p0){
          if (sum(pred.cv1$proba > p0)){
            paths1 <- pred.cv1$paths[pred.cv1$proba > p0][pred.cv1$beta.df[[folds.ncol[which(p0.grid == p0),fold1]]] > 0]
          }else{
            paths1 <- list()
          }
          if (sum(pred.cv2$proba > p0)){
            paths2 <- pred.cv2$paths[pred.cv2$proba > p0][pred.cv2$beta.df[[folds.ncol[which(p0.grid == p0),fold2]]] > 0]
          }else{
            paths2 <- list()
          }
          len <- (length(paths1) + length(paths2))/2
          if (len > 0){
            c(length(intersect(paths1, paths2))/len, len)
          }else{
            c(1, 0)
          }
        })
      })
    })
    num.rule.df <- lapply(1:(nfold - 1), function(fold1){
      sapply(1:(nfold-fold1), function(fold2){
        stab.df[[fold1]][[fold2]][2,]
      })
    })
    num.rule.df <- do.call('cbind', num.rule.df)
    num.rules <- apply(num.rule.df, 1, mean)
    stab.df <- lapply(1:(nfold - 1), function(fold1){
      sapply(1:(nfold-fold1), function(fold2){
        stab.df[[fold1]][[fold2]][1,]
      })
    })
    stab.df <- do.call('cbind', stab.df)
    stab.grid <- apply(stab.df, 1, mean)

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

  ind.min <- ind.1 - 1 + which.min(error.mean[ind.1:ind.max])
  ind.pred <- max(ind.1, min(which(error.mean <= error.mean[ind.min] + 2*error.sd[ind.min])))
  p0.stab.cv <- sapply(error.grids, function(grid.cv){
    criterion <- (grid.cv$stab.grid[ind.1:ind.max] - 0.90)^2 + (0.0 - grid.cv$error.grid[ind.1:ind.max])^2
    p0.grid[ind.1:ind.max][which.min(criterion)]
  })
  p0.stab <- median(p0.stab.cv)
  
  return(list(p0.pred = p0.grid[ind.pred], p0.stab = p0.stab, error.grid.p0 = error.grid.p0, type = type))
  
}

#'
#' Plot SIRUS cross-validation path: error and stability versus the number of rules when \code{p0} varies.
#'
#' Error is 1-AUC for classification and the unexplained variance for regression.
#' Stability is the average proportion of rules shared by two SIRUS models fit on two distinct folds of the cross-validation.
#'
#' @title Plot SIRUS cross-validation path.
#'
#' @param sirus.cv.grid Cross-validation results returned by \code{\link{sirus.cv}}.
#' @param p0.criterion Criterion to pick the optimal \code{p0} displayed in the plots: if 'pred' then \code{p0.pred} is used for a minimal error, if 'stab' then \code{p0.stab} is used for a tradeoff error/stability. Default is 'pred' for classification and 'stab' for regression.
#' @param num.rule.max Upper limit on the number of rules for the x-axis. Default is 25.
#'
#' @return Plots of cross-validation results.
#'   \item{\code{error}}{plot of error vs number of rules (ggplot2 object).}
#'   \item{\code{stability}}{plot of stability vs number of rules (ggplot2 object).}
#' @export
#'
#' @examples
#' ## load SIRUS
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
sirus.plot.cv <- function(sirus.cv.grid, p0.criterion = NULL, num.rule.max = 25){
  
  # check arguments
  # criterion for optimal p0
  if (is.null(p0.criterion)){
    if (sirus.cv.grid$type == 'reg'){p0.criterion <- 'stab'}else{p0.criterion <- 'pred'}
  }else{
    if (!p0.criterion %in% c('pred', 'stab')){
      stop('Invalid p0 criterion. Its value has to be pred or stab.')
    }
  }
  # max num of rules
  num.rule.valid <- is.numeric(num.rule.max) & num.rule.max > 0
  if (!num.rule.valid){
    stop('Invalid maximum number of rules. Its value should be a positive integer.')
  }

  # filter cv grid & performance metrics
  error.grid.p0 <- sirus.cv.grid$error.grid.p0
  num.rule.max <- min(num.rule.max, max(error.grid.p0$num.rules.mean))
  grid.index <- sapply(1:num.rule.max, function(ind){which.min(abs(ind - sirus.cv.grid$error.grid.p0$num.rules.mean))[1]})
  if (p0.criterion == 'stab'){
    p0 <- sirus.cv.grid$p0.stab
  }else{
    p0 <- sirus.cv.grid$p0.pred
  }
  ind.p0 <- which.min(abs(sirus.cv.grid$error.grid.p0$p0.grid - p0))
  grid.cv <- sirus.cv.grid$error.grid.p0[ind.p0,]
  num.rules <- grid.cv$num.rules.mean
  error <- grid.cv$error.mean
  stab <- grid.cv$stab.mean
  error.grid.p0 <- sirus.cv.grid$error.grid.p0[sort(unique(c(grid.index,ind.p0))),]
  
  # declare variables
  num.rules.mean <- NULL
  error.mean <- NULL
  error.sd <- NULL
  stab.mean <- NULL
  stab.sd <- NULL

  # plot error vs number of rules
  label.error <- if (sirus.cv.grid$type == 'reg'){'Unexplained variance'}else{'1-AUC'}
  tag.y <- (min(error.grid.p0$error.mean) + max(error.grid.p0$error.mean))/2
  plot.error <- ggplot(error.grid.p0, aes(x = num.rules.mean, y = error.mean)) +
    geom_line(size = 0.8) + geom_point(size=1) +
    geom_errorbar(aes(ymin=error.mean - 2*error.sd,
                      ymax=error.mean + 2*error.sd),  width=0.7, size=0.5)+
    geom_hline(yintercept = error,
               linetype = 'dashed', color = 'blue', size = 0.7) +
    geom_vline(xintercept = num.rules,
               linetype = 'dashed', color = 'blue', size = 0.7) +
    geom_text(aes(num.rules, tag.y, label = 'Optimal p0'), angle = '90',
              vjust=1.5, color='blue', size = 7) +
    geom_text(aes(num.rule.max - 8, error, label = paste0('SIRUS Error: ', round(error,2))),
              vjust=-1, color='blue', size = 7) +
    xlab('Number of rules') +
    ylab(label.error) +
    theme_classic() +
    theme(axis.line.x = element_line(colour = 'black'),
          axis.line.y = element_line(colour = 'black'),
          text = element_text(size=20),
          plot.title = element_text(hjust = 0.5, size=18, face="italic"))

  # plot stability vs number of rules
  stab.extremum <- c(min(error.grid.p0$stab.mean), max(error.grid.p0$stab.mean))
  tag.p0 <- (stab.extremum[which.max(abs(stab.extremum - stab))] + stab)/2
  plot.stab <- ggplot(error.grid.p0, aes(x = num.rules.mean, y = stab.mean)) +
    geom_line(size = 0.8) + geom_point() +
    geom_errorbar(aes(ymin=stab.mean - stab.sd,
                      ymax=stab.mean + stab.sd),  width=0.7, size=0.5)+
    geom_hline(yintercept = stab,
               linetype = 'dashed', color = 'blue', size = 0.7) +
    geom_vline(xintercept = num.rules,
               linetype = 'dashed', color = 'blue', size = 0.7) +
    geom_text(aes(num.rules, tag.p0, label = 'Optimal p0'), angle = '90',
              vjust=1.5, color='blue', size = 7) +
    geom_text(aes(num.rules - 5, stab, label = paste0('SIRUS Stability: ', round(stab,2))),
              vjust=-1, color='blue', size = 7) +
    xlab('Number of rules') +
    ylab('Stability') +
    theme_classic() + ylim(c(0.9*min(error.grid.p0$stab.mean), 1.1*max(error.grid.p0$stab.mean))) +
    theme(axis.line.x = element_line(colour = 'black'),
        axis.line.y = element_line(colour = 'black'),
        text = element_text(size=20),
        plot.title = element_text(hjust = 0.5, size=18, face="italic"))

  return(list(error = plot.error, stability = plot.stab))

}
