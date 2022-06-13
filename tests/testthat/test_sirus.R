library(sirus)
context("sirus")

test_that("sirus.fit returns the defined number of rules.", {
  data <- iris
  y <- rep(0, nrow(data))
  y[data$Species == 'setosa'] = 1
  data$Species <- NULL
  num.rule <- 8
  sirus.m <- sirus.fit(data, y, num.rule = num.rule, p0 = NULL, num.rule.max = 20, q = 10,
                       num.trees.step = 1000, alpha = 0.05, mtry = NULL, num.trees = NULL,
                       num.threads = NULL, replace = TRUE, verbose = FALSE, seed = NULL)
  model.size <- length(sirus.m$rules)
  expect_equal(model.size, num.rule)
})

test_that("sirus.print returns the defined number of rules.", {
  data <- iris
  y <- rep(0, nrow(data))
  y[data$Species == 'setosa'] = 1
  data$Species <- NULL
  num.rule <- 8
  sirus.m <- sirus.fit(data, y, num.rule = num.rule, p0 = NULL, num.rule.max = 20, q = 10,
                       num.trees.step = 1000, alpha = 0.05, mtry = NULL, num.trees = NULL,
                       num.threads = NULL, replace = TRUE, verbose = FALSE, seed = NULL)
  m.print <- sirus.print(sirus.m)
  model.size <- length(m.print) - 1
  expect_equal(model.size, num.rule)
})

test_that("sirus.predict returns valid predictions.", {
  data <- iris
  y <- rep(0, nrow(data))
  y[data$Species == 'setosa'] = 1
  data$Species <- NULL
  num.rule <- 8
  sirus.m <- sirus.fit(data, y, num.rule = num.rule, p0 = NULL, num.rule.max = 20, q = 10,
                       num.trees.step = 1000, alpha = 0.05, mtry = NULL, num.trees = NULL,
                       num.threads = NULL, replace = TRUE, verbose = FALSE, seed = NULL)
  predictions <- sirus.predict(sirus.m, data)
  is.valid <- is.numeric(predictions) & all(predictions >= 0) & all(predictions <= 1)
  expect_true(is.valid)
})

test_that("sirus.cv returns valid p0", {
  data <- iris
  y <- rep(0, nrow(data))
  y[data$Species == 'setosa'] = 1
  data$Species <- NULL
  sirus.cv.grid <- sirus.cv(data, y, nfold = 3, ncv = 2, num.rule.max = 25, q = 10,
                       num.trees.step = 1000, alpha = 0.05, mtry = NULL, num.trees = 100,
                       num.threads = NULL, replace = TRUE, sample.fraction = NULL,
                       verbose = FALSE, seed = NULL)
  p0 <- sirus.cv.grid$p0.stab
  is.valid <- is.numeric(p0) & (p0 > 0) & (p0 < 1)
  expect_true(is.valid)
})

