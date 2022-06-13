
# Merge paths and proba of two ranger runs
rangerMerge <- function(numTrees, paths1, proba1, paths2, proba2){
  res <- rangerMergeCpp(numTrees, paths1, proba1, paths2, proba2)
  return(res)
}

# compute stability metric
stabilityMetric <- function(numTrees, proba){
  res <- stabilityMetricCpp(numTrees, proba)
  return(res)
}

# bin variable with empirical q-quantiles
get.bins <- function(X, y, q, discrete.limit){
  bins <- list()
  if (is.numeric(X)){
    values <- unique(X)
    if (length(values) > discrete.limit){
      # continuous variable
      bins$type <- 'continuous'
      bins$cut.values <- quantile(X, probs = seq(0,1,1/q))
      names(bins$cut.values) <- NULL
    }else{
      # discrete variable
      bins$type <- 'discrete'
      bins$cut.values <- sort(values)
    }
  }
  if (is.factor(X)){
    # categorical variable
    bins$type <- 'categorical'
    levels.order <- order(sapply(levels(X), function(level){
      mean(y[X == level])
    }))
    bins$levels <- levels(X)[levels.order]
  }
  return(bins)
}
binarize.X <- function(X, bins, q){
  # continuous and discrete variables
  if (bins$type %in% c('continuous', 'discrete')){
    breaks <- bins$cut.values
    cut.on.max <- (diff(breaks[length(breaks) + -1:0]) == 0 || bins$type == 'discrete') # cuts on maximum values are possible only for discrete variables
    breaks <- unique(breaks)
    if (length(breaks) >= 2){
      X.bin <- cut(as.numeric(X), breaks = breaks, labels = FALSE, right = FALSE)
      if (cut.on.max){
        X.bin[is.na(X.bin)] <- length(breaks)
      }else{
        X.bin[is.na(X.bin)] <- length(breaks) - 1
      }
    }else{
      # constant variable
      X.bin <- rep(1, length(X))
    }
  }
  # categorical variable
  if (bins$type == 'categorical'){
    X.bin <- factor(X, bins$levels) # order levels
    X.bin <- as.numeric(X.bin) # transform categorical in ordered variables
  }
  return(X.bin)
}

# run ranger iteratively
ranger.stab <- function(data.bin.y, num.trees.step = 1000, alpha = 0.05, mtry = NULL, max.depth = 2, 
                        num.trees = NULL, num.threads = NULL, replace = TRUE,
                        sample.fraction = ifelse(replace, 1, 0.632), verbose = TRUE, seed = NULL){
  
  if (is.null(num.trees)){
    # automatic set of num.trees with stability stopping criterion
    stab.metric <- 0
    num.trees <- 0
    # Build trees until 1 - alpha (95%) stability
    while (stab.metric < (1 - alpha)){
      forest <- ranger(y ~., data = data.bin.y, num.trees = num.trees.step, mtry = mtry, max.depth = max.depth,
                       importance = 'none', oob.error = FALSE, write.forest = FALSE, replace = replace,
                       sample.fraction = sample.fraction, num.threads = num.threads, seed = seed)
      paths.temp <- forest$paths
      proba.temp <- forest$paths.proba
      num.trees <- num.trees + num.trees.step
      if (num.trees == num.trees.step){
        paths <- paths.temp
        proba <- proba.temp
        stab.metric <- stabilityMetric(num.trees, proba/num.trees)
      }else{
        # merge paths and probas
        paths.merge <- rangerMerge(num.trees, paths, proba, paths.temp, proba.temp)
        paths <- paths.merge[[1]]
        proba <- paths.merge[[2]]
        stab.metric <- paths.merge[[3]]
      }
      if (verbose){
        print(paste('Number of trees: ', num.trees, ' - Stability ', 100*round(stab.metric, 4), ' %.'))
      }
    }
  }else{
    # User provided value for num.trees
    forest <- ranger(y ~., data = data.bin.y, num.trees = num.trees, mtry = mtry, max.depth = max.depth,
                     importance = 'none', oob.error = FALSE, write.forest = FALSE, replace = replace,
                     sample.fraction = sample.fraction, num.threads = num.threads, seed = seed)
    paths <- forest$paths
    proba <- forest$paths.proba
  }
  
  # format paths and proba
  proba <- proba/num.trees
  index.proba <- order(-proba)
  paths <- paths[index.proba]
  proba <- proba[index.proba]
  
  return(list(paths = paths, proba = proba, num.trees = num.trees))
  
}

# path post-treatment for d = 1 or 2 (exact and deterministic algorithm)
paths.filter.2 <- function(paths, proba, num.rule){
  
  paths.ftr <- list()
  proba.ftr <- c()
  split.gen <- list()
  ind.max <- length(paths)
  ind <- 1
  num.rule.temp <- 0
  
  while (num.rule.temp < (num.rule) & ind <= ind.max){
    
    path.ind <- paths[[ind]]
    
    # Remove empty split (variable 0)
    split.var <- sapply(path.ind, function(split){
      split[1]
    })
    if (0 %in% split.var){
      path.ind <- path.ind[which(split.var != 0)]
    }
    
    # format rule with 2 cuts on same variable
    ## a mettre ds le C++
    if (length(path.ind) == 2){
      if (path.ind[[1]][1] == path.ind[[2]][1] & path.ind[[1]][3] == path.ind[[2]][3]){
        if (path.ind[[1]][2] > path.ind[[2]][2]){
          if (path.ind[[1]][3] == 1){
            path.ind <- path.ind[1]
          }else{
            path.ind <- path.ind[2]
          }
        }else{
          if (path.ind[[1]][3] == 1){
            path.ind <- path.ind[2]
          }else{
            path.ind <- path.ind[1]
          }
        }
        paths[[ind]] <- path.ind
      }
    }

    split.ind <- lapply(path.ind, function(split){
      split[1:2]
    })
    d <- length(path.ind)
    
    if (!list(split.ind) %in% split.gen){
      
      ### add path
      paths.ftr <- c(paths.ftr, list(path.ind))
      proba.ftr <- c(proba.ftr, proba[ind])
      num.rule.temp <- length(paths.ftr)
      
      ### add generated interaction
      if (d <= 2){
        # 1-split path
        # TO DO : add path 2 splits same var
        if (d == 1){
          split.gen.temp <- list(split.ind)
          split.gen <- c(split.gen, split.gen.temp)
        }
        # 2 split-path
        if (d == 2){
          
          # get index of rules involving any similar constraint
          bool.ind <- lapply(paths.ftr, function(path){
            if (length(path) <= 2){
              bool <- sapply(path, function(x){
                any(sapply(path.ind, function(y){
                  all(y[1:2] == x[1:2])
                }))
              })
              return(c(all(bool), any(bool)))
            }else{
              return(c(FALSE, FALSE))
            }
          })
          bool.all <- sapply(bool.ind, function(x){x[1]})
          bool.any <- sapply(bool.ind, function(x){x[2]})
          bool.mixed <- !bool.all & bool.any
          num.rule.all <- sum(bool.all)
          num.rule.any <- sum(bool.any)
          
          if (num.rule.all >= 2){
            split.gen <- c(split.gen, list(split.ind))
          }
          
          # combine path with paths.ftr
          split.gen.temp <- lapply(paths.ftr[bool.mixed], function(x){
            split.diff <- setdiff(c(x, path.ind), intersect(x, path.ind))
            split1 <- list(list(split.diff[[1]][1:2]))
            if (all(split.diff[[1]][1:2] == split.diff[[2]][1:2]) & !(split1 %in% split.gen)){
              return(split1)
            }
          })
          # specFific case of two splits on the same direction
          if (split.ind[[1]][1] == split.ind[[2]][1]){
            bool.double <- sapply(split.gen, function(split){
              all(sapply(split, function(x){list(x) %in% split.ind})) & length(split) == 1
            })
            if (length(bool.double) > 0){
              split.gen.temp2 <- lapply(split.gen[bool.double], function(split){
                split.diff <- setdiff(split.ind, split)
                if (length(split.diff) > 0){
                  split.diff
                }
              })
              split.gen.temp <- c(split.gen.temp, split.gen.temp2)
            }
          }
        }
        
        split.gen.temp <- Filter(Negate(is.null), split.gen.temp)
        if (length(split.gen.temp) > 0){
          split.gen.1 <- Filter(function(x){length(x)==1}, split.gen)
          split.gen.temp <- c(split.gen.temp, unlist(lapply(split.gen.temp, function(split){
            lapply(split.gen.1, function(split1){
              if (length(split) == 1 & split[[1]][1] == split1[[1]][1] & split[[1]][2] != split1[[1]][2]){
                if (split[[1]][2] > split1[[1]][2]){
                  c(split1, split)
                }else{
                  c(split, split1)
                }
              }
            })
          }), recursive = F))
          split.gen.temp <- setdiff(split.gen.temp, split.gen)
          split.gen <- c(split.gen, split.gen.temp)
        }
      }
    }
    ind <- ind + 1
    
  }
  
  return(list(paths = paths.ftr, proba = proba.ftr))
    
}
# path post-treatment  for any d (stochastic heuristic)
paths.filter.d <- function(paths, proba, num.rule, data.bin){
  
  # parameters
  nvar <- ncol(data.bin)
  nrule.max <- length(paths)
  nrule <- 1
  paths.ftr <- list()
  proba.ftr <- list()
  ind.batch <- 0
  nsample <- 10*num.rule
  
  # generate independent data
  values <- lapply(1:nvar, function(j){unique(data.bin[,j])})
  data.indep <- sapply(values, function(x){sample(x, size = nsample, replace = T)})
  
  while (nrule < num.rule & ind.batch < nrule.max){
    
    # select a batch of rules
    nrule.batch <- min(nrule.max - ind.batch, num.rule)
    paths.k <- paths[ind.batch + 1:nrule.batch]
    proba.k <- proba[ind.batch + 1:nrule.batch]
    
    # get data.rule
    data.rule.b <- sapply(paths.k, get.rule.support.from.bin, data.bin = data.indep)
    if (ind.batch == 0){
      ind.sel <- c(1)
      data.rule <- cbind(rep(1, nsample), data.rule.b)
    }else{
      ind.sel <- 1:ncol(data.rule)
      data.rule <- cbind(data.rule, data.rule.b)
    }
    
    # get index of rules to remove
    ind.path <- list()
    for (ind in (nrule + 1):ncol(data.rule)){
      ind.sel <- c(ind.sel, ind)
      data.rule.temp <- data.rule[, ind.sel, drop = F]
      rk.diff <- length(ind.sel) - rankMatrix(data.rule.temp)[1]
      rk.diff
      if (rk.diff > 0){ind.sel <- head(ind.sel, -1)}
    }
    
    data.rule <- data.rule[, ind.sel, drop = F]
    if (ind.batch == 0){ind.path <- ind.sel - nrule}else{ind.path <- ind.sel - nrule - 1}
    ind.path <- ind.path[ind.path > 0]
    paths.ftr <- c(paths.ftr, paths.k[unlist(ind.path)])
    proba.ftr <- c(proba.ftr, proba.k[unlist(ind.path)])
    nrule <- length(paths.ftr)
    ind.batch <- ind.batch + nrule.batch
    
  }
  
  paths.ftr <- paths.ftr[1:min(num.rule, nrule)]
  proba.ftr <- unlist(proba.ftr[1:min(num.rule, nrule)])
  
  return(list(paths = paths.ftr, proba = proba.ftr))
  
}

# format path
format.path <- function(path, bins.list){
  lapply(path, function(split){
    bins <- bins.list[[split[1]]]
    if (bins$type %in% c('continuous', 'discrete')){
      # continuous and discrete variables
      breaks <- bins$cut.values
      split[2] <- min(which(breaks == unique(breaks)[ceiling(split[2])])) - 1
      split[3] <- if(split[3] == 0){'L'}else{'R'}
      split[4] <- breaks[as.numeric(split[2]) + 1]
    }
    if (bins$type == 'categorical'){
      # categorical variables
      if (split[3] == 0){
        split <- c(split[1], '=', sort(bins$levels[1:floor(as.numeric(split[2]))]))
      }else{
        split <- c(split[1], '=', sort(bins$levels[ceiling(as.numeric(split[2])):length(bins$levels)]))
      }
    }
    return(split)
  })
}

# recover rule from path
get.rule <- function(path, bins.list, data.names){
  lapply(path, function(split){
    var.name <- data.names[as.numeric(split[1])]
    split.type <- bins.list[[as.numeric(split[1])]]$type
    if (split.type %in% c('continuous', 'discrete')){
      if(split[3]=='L'){sign <- '<'}
      if(split[3]=='R'){sign <- '>='}
      split.value <- split[4]
      rule <- c(var.name, sign, split.value)
    }
    if (split.type %in% c('categorical')){
      rule <- split
      rule[1] <- var.name
    }
    rule
  })
}

# get rule outputs
get.rule.support.from.bin <- function(path, data.bin){
  splits <- sapply(path, function(split){
    X <- data.bin[, split[1]]
    if (split[3] == 0){
      X <- X < split[2]
    }
    if (split[3] == 1){
      X <- X >= split[2]
    }
    X
  })
  apply(splits, 1, prod)
}
get.rule.support <- function(data, rules){  
  as.data.frame(lapply(rules, function(rule){
    Z <- sapply(rule, function(split){
      X <- data[,split[1]]
      if (is.numeric(X)){
        if (split[2] == '<'){
          Z <- X < as.numeric(split[3])
        }
        if (split[2] == '>='){
          Z <- X >= as.numeric(split[3])
        }
      }
      if (is.factor(X)){
        Z <- X %in% split[3:length(split)]
      }
      Z
    })
    if (!is.matrix(Z)){
      all(Z)
    }else{
      apply(Z, 1, all)
    }
  }))
}
get.rule.outputs <- function(data.rule.supp, y){
  lapply(1:ncol(data.rule.supp), function(j){
    Z <- data.rule.supp[, j]
    outputs <- c(mean(y[Z == 1]), mean(y[Z == 0]))
    supp.size <- c(sum(Z), sum(1 - Z))
    list(outputs = outputs, supp.size = supp.size)
  })
}

# transform data
get.data.rule <- function(data, rules, rules.out){
  rules.bool <- get.rule.support(data, rules)
  ndata <- nrow(data)
  data.rule <- matrix(sapply(1:length(rules), function(ind) {
    if (ndata > 1){
      rule.out <- rep(rules.out[[ind]]$outputs[2], ndata)
      rule.out[rules.bool[,ind]] <- rules.out[[ind]]$outputs[1]
    }else{
      rule.out <- if(rules.bool[,ind]){rules.out[[ind]]$outputs[1]}else{rules.out[[ind]]$outputs[2]}
    }
    rule.out
  }), nrow = ndata)
  return(data.rule)
}

# check data
data.check <- function(data){
  
  # check data type
  is.df <- is.data.frame(data)
  if (!is.df){
    stop('Invalid data. data is not a dataframe.')
  }else{
    if (!all(dim(data) > 0)){
      stop('Invalid data. data is empty.')
    }
  }
  # check data is either numeric or factor
  is.type <- sapply(1:ncol(data), function(ind){
    is.numeric(data[,ind]) || is.factor(data[,ind])
  })
  if (!all(is.type)){
    offending_columns <- colnames(data)[!is.type]
    stop("Invalid data. Variable types have to be numeric or factor. Not the case for columns: ",
         paste0(offending_columns, collapse = ", "), ". ", call. = FALSE)
  }
  # check missing values
  if (any(is.na(data))) {
    offending_columns <- colnames(data)[colSums(is.na(data)) > 0]
    stop("Invalid data. Missing data in columns: ",
         paste0(offending_columns, collapse = ", "), ".", call. = FALSE)
  }
  
}

# check data and y
data.y.check <- function(data, y){
  
  # check data
  data.check(data)
  
  # check output
  y.valid <- is.vector(y) & is.numeric(y)
  if (!y.valid){
    stop("Invalid y. The output y must be a numeric vector.")
  }
  # check dimension consistency
  dim.valid <- nrow(data) == length(y)
  if (!dim.valid){
    stop('Invalid data and y. data and y should have compatible dimensions.')
  }
  # check sample size
  size.valid <- nrow(data) > 5
  if (!size.valid){
    warning('Small sample size. At least 6 data points are required to return a non-empty rule list.')
  }
  
}

# check sirus parameters
sirus.param.check <- function(data, num.rule.max, q, num.trees.step, alpha, mtry){
  
  # check num.rule.max
  num.rule.max.valid <- is.numeric(num.rule.max)
  if (num.rule.max.valid){num.rule.max.valid <- (round(num.rule.max) == num.rule.max) & num.rule.max > 0}
  if (!num.rule.max.valid){
    stop("Invalid num.rule.max. Maximum number of rules must be a positive integer.")
  }else{
    if (num.rule.max > 100){
      warning('Warning num.rule.max: SIRUS is designed to output short list of rules (typically < 100 rules).')
    }
  }
  # check q
  q.valid <- is.numeric(q)
  if (q.valid){q.valid <- (round(q) == q) & q >= 2}
  if (!q.valid){
    stop("Invalid q. Number of quantiles must be an integer greater or than to 2.")
  }else{
    if (q > nrow(data)){
      warning('Warning q: Number of quantiles should be much smaller than the sample size.')
    }
  }
  # alpha
  alpha.valid <- all(c(is.numeric(alpha), alpha > 0, alpha < 1))
  if (!alpha.valid){
    stop('Invalid alpha. alpha must be a numeric value between 0 and 1.')
  }
  # mtry
  if (!is.null(mtry)){
    mtry.valid <- is.numeric(mtry)
    if (mtry.valid){mtry.valid <- (round(mtry) == mtry) & mtry > 0 & mtry <= ncol(data)}
    if (!mtry.valid){
      stop("Invalid mtry. mtry must be a positive integer smaller than the number of variables.")
    }
  }
  # num.tree.step
  num.trees.step.valid <- is.numeric(num.trees.step)
  if (num.trees.step.valid){num.trees.step.valid <- (round(num.trees.step) == num.trees.step) & num.trees.step >= 10}
  if (!num.trees.step.valid){
    stop("Invalid num.trees.step. The step of the number of trees must be an integer greater than 10.")
  }
  
}

# check sirus model
sirus.model.check <- function(sirus.m){
  
  sirus.valid <- FALSE
  type.valid <- is.list(sirus.m)
  if (type.valid){
    names.valid <- all(names(sirus.m) == c('rules', 'rules.out', 'proba', 'paths', 
                      'rule.weights', 'rule.glm', 'type', 'num.trees', 'data.names', 'mean', 'bins'))
    if (names.valid){
      rules.valid <- all(sapply(sirus.m$rules, function(rule){
        is.list(rule) & 
        all(sapply(rule, function(split){
          if (split[2] == '='){length(split) >= 3}else{length(split) == 3}
        }))
      }))
      rules.out.valid <- all(sapply(sirus.m$rules.out, function(rule.out){
        is.list(rule.out) &
        all(names(rule.out) == c('outputs', 'supp.size')) &
        is.numeric(unlist(rule.out))
      }))
      proba.valid <- is.numeric(sirus.m$proba) & all(sirus.m$proba >= 0) & all(sirus.m$proba <= 1)
      mean.valid <- is.numeric(sirus.m$mean)
      weights.valid <- is.numeric(sirus.m$rule.weights)
      sirus.valid <- type.valid & names.valid & rules.valid & rules.out.valid &
                     proba.valid & mean.valid
    }
  }
  return(sirus.valid)
  
}
