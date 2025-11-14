# Data Cleaning, Diagnosis, Imputation, Transformation,
# Dimensionality Reduction, Frequent Patterns, Sequence Mining,
# Statistical Tests, and Clustering — Reusable R Script

# -----------------------------
# 0. Setup: install/load packages
# -----------------------------
.pkgs <- c(
  'tidyverse', 'data.table', 'naniar', 'VIM', 'janitor', 'lubridate',
  'scales', 'caret', 'factoextra', 'cluster', 'dbscan', 'Rtsne', 'uwot',
  'arules', 'arulesViz', 'arulesSequences', 'broom'
)
# install missing packages (uncomment to install)
.inst <- .pkgs[!(.pkgs %in% installed.packages()[, 'Package'])]
if(length(.inst)) message('To install missing packages, run: install.packages(', paste0("'", .inst, "'", collapse = ", "), ")")

lapply(.pkgs, function(p) suppressPackageStartupMessages(require(p, character.only = TRUE)))

# -----------------------------
# 1. Diagnosis functions
# -----------------------------

# Summarize dataset structure and basic issues
diagnose_data <- function(df) {
  stopifnot(is.data.frame(df))
  n <- nrow(df); p <- ncol(df)
  cat(sprintf("Rows: %d\nColumns: %d\n", n, p))
  cat('\n-- Column types and missingness --\n')
  types <- sapply(df, class)
  miss <- sapply(df, function(x) sum(is.na(x)))
  pct_miss <- round(miss / n * 100, 2)
  summary_df <- data.frame(column = names(df), type = types, missing = miss, pct_missing = pct_miss)
  print(as_tibble(summary_df))
  cat('\n-- Duplicate rows --\n')
  cat(sprintf('Duplicate rows: %d\n', sum(duplicated(df))))
  cat('\n-- Basic numeric summary (first 10 numeric cols) --\n')
  num <- df %>% select(where(is.numeric))
  if(ncol(num)>0) print(summary(num)) else cat('No numeric columns\n')
  invisible(summary_df)
}

# Visual missingness using naniar
plot_missingness <- function(df) {
  if(!requireNamespace('naniar', quietly = TRUE)) stop('Please install naniar to use plot_missingness')
  naniar::gg_miss_upset(df)
}

# Detect special placeholder values often used to indicate missing or special
detect_special_values <- function(df, specials = c(-99, -9999, 9999, 'NA', 'N/A', ''), top_n = 10) {
  res <- list()
  for(col in names(df)){
    vals <- df[[col]]
    # compare after coercion to character
    svals <- unique(as.character(vals))
    found <- intersect(as.character(specials), svals)
    if(length(found)) res[[col]] <- found
  }
  res
}

# Detect outliers: IQR and Z-score methods
detect_outliers_iqr <- function(x, na.rm = TRUE, coef = 1.5) {
  if(!is.numeric(x)) stop('x must be numeric')
  x2 <- x[!is.na(x)]
  Q1 <- quantile(x2, 0.25)
  Q3 <- quantile(x2, 0.75)
  IQR <- Q3 - Q1
  lower <- Q1 - coef * IQR
  upper <- Q3 + coef * IQR
  which(x < lower | x > upper)
}

detect_outliers_zscore <- function(x, threshold = 3, na.rm = TRUE) {
  x <- as.numeric(x)
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  which(abs((x - m)/s) > threshold)
}

# Inconsistency detection: mixed types or unexpected levels
detect_inconsistencies <- function(df) {
  problems <- list()
  for(col in names(df)){
    x <- df[[col]]
    # mixed types
    types <- unique(sapply(x, function(v) class(v)[1]))
    if(length(types) > 1) problems[[col]] <- list(type_issue = types)
    # factors/characters with many unique values could be inconsistent
    if(is.character(x) || is.factor(x)){
      uniq <- unique(x)
      if(length(uniq) < length(x) && any(grepl('\\s+$|^\\s+|\\s{2,}', uniq))){
        problems[[col]] <- c(problems[[col]], list(whitespace_levels = uniq[grepl('\\s', uniq)]))
      }
    }
  }
  problems
}

# Localization checks: detect date-like columns, decimal separators, and thousand separators
detect_localization_issues <- function(df, sample_n = 100) {
  res <- list()
  for(col in names(df)){
    x <- na.omit(as.character(df[[col]]))
    if(length(x)==0) next
    s <- head(x, sample_n)
    # date detection
    dt_parsed <- suppressWarnings(lubridate::parse_date_time(s, orders = c('ymd', 'dmy', 'mdy','Ymd HMS', 'dmY'))) 
    pct_dates <- mean(!is.na(dt_parsed))
    # decimal separators: presence of comma and dot
    pct_comma <- mean(grepl(',', s) & !grepl('\s', s))
    pct_dot <- mean(grepl('\.', s))
    if(pct_dates > 0.3) res[[col]] <- c(res[[col]], date_like = pct_dates)
    if(pct_comma > 0.3) res[[col]] <- c(res[[col]], comma_decimal = pct_comma)
    if(pct_dot > 0.3) res[[col]] <- c(res[[col]], dot_decimal = pct_dot)
  }
  res
}

# -----------------------------
# 2. Transformations and pipelines
# -----------------------------

# Common transformations pipeline using caret preProcess
transform_pipeline <- function(df, numeric_cols = NULL, center = TRUE, scale = TRUE, BoxCox = FALSE) {
  if(is.null(numeric_cols)) numeric_cols <- names(df)[sapply(df, is.numeric)]
  pre <- caret::preProcess(df[, numeric_cols, drop = FALSE], method = c(if(center) 'center' else NULL, if(scale) 'scale' else NULL, if(BoxCox) 'BoxCox' else NULL))
  df2 <- df
  df2[, numeric_cols] <- predict(pre, df[, numeric_cols, drop = FALSE])
  list(data = df2, preproc = pre)
}

# Log-safe transform
safe_log <- function(x, offset = 1e-6){
  x2 <- ifelse(x <= 0 | is.na(x), NA, x)
  log(x2 + offset)
}

# -----------------------------
# 3. Deductive correction and deterministic imputation
# -----------------------------

# Apply user-specified deterministic rules. rules is a list of functions: function(df) df
apply_deductive_corrections <- function(df, rules = list()){
  for(i in seq_along(rules)){
    rule <- rules[[i]]
    if(!is.function(rule)) next
    df <- rule(df)
  }
  df
}

# Deterministic imputation helpers
# - mode imputation for categorical
# - median/mean for numeric
# - last observation carried forward (for time series)

impute_mode <- function(x){
  ux <- na.omit(unique(x))
  if(length(ux) == 0) return(x)
  tab <- sort(table(x), decreasing = TRUE)
  modev <- names(tab)[1]
  x[is.na(x)] <- modev
  x
}

impute_numeric <- function(x, method = c('median', 'mean'), by = NULL){
  method <- match.arg(method)
  if(!is.null(by)){
    # by is grouping vector
    df <- data.frame(x = x, by = by)
    df <- df %>% group_by(by) %>% mutate(x = ifelse(is.na(x), if(method=='median') median(x, na.rm=TRUE) else mean(x, na.rm=TRUE), x)) %>% ungroup()
    return(df$x)
  }
  if(method=='median') x[is.na(x)] <- median(x, na.rm = TRUE) else x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
}

locf <- function(x) {
  if(!requireNamespace('zoo', quietly = TRUE)) stop('Install zoo for locf/backfill')
  zoo::na.locf(x, na.rm = FALSE)
}

# -----------------------------
# 4. Dimensionality reduction
# -----------------------------

run_pca <- function(df, numeric_cols = NULL, scale. = TRUE, center = TRUE, ncomp = NULL){
  if(is.null(numeric_cols)) numeric_cols <- names(df)[sapply(df, is.numeric)]
  mat <- na.omit(df[, numeric_cols, drop = FALSE])
  pca <- prcomp(mat, center = center, scale. = scale.)
  if(is.null(ncomp)) ncomp <- min(ncol(mat), 5)
  list(pca = pca, scores = pca$x[,1:ncomp, drop = FALSE])
}

run_tsne <- function(df, numeric_cols = NULL, perplexity = 30, theta = 0.5, max_iter = 1000){
  if(is.null(numeric_cols)) numeric_cols <- names(df)[sapply(df, is.numeric)]
  mat <- na.omit(as.matrix(df[, numeric_cols, drop = FALSE]))
  Rtsne::Rtsne(mat, perplexity = perplexity, theta = theta, max_iter = max_iter)
}

run_umap <- function(df, numeric_cols = NULL, n_neighbors = 15, n_components = 2){
  if(is.null(numeric_cols)) numeric_cols <- names(df)[sapply(df, is.numeric)]
  mat <- na.omit(as.matrix(df[, numeric_cols, drop = FALSE]))
  uwot::umap(mat, n_neighbors = n_neighbors, n_components = n_components)
}

# -----------------------------
# 5. Frequent patterns & sequence mining
# -----------------------------

# Frequent itemsets with arules (transactions input expected)
frequent_itemsets <- function(transactions, support = 0.01, maxlen = 10){
  if(is.data.frame(transactions)) transactions <- as(split(transactions$item, transactions$transactionID), 'transactions')
  apriori(transactions, parameter = list(supp = support, target = 'frequent', maxlen = maxlen))
}

association_rules <- function(transactions, support = 0.01, confidence = 0.5){
  apriori(transactions, parameter = list(supp = support, conf = confidence, target = 'rules'))
}

# Sequence mining: arulesSequences expects a special sparse data format
sequence_mining <- function(df_sequences, support = 0.01, maxsize = 10){
  # df_sequences should be a data.frame with columns: sequenceID, eventID, item
  if(!('sequenceID' %in% names(df_sequences) && 'eventID' %in% names(df_sequences) && 'item' %in% names(df_sequences))) stop('df_sequences must have sequenceID, eventID, and item columns')
  seqdata <- as(df_sequences, 'transactions') # quick attempt; for arulesSequences, one usually uses cspade on transactions-like sequences
  # If arulesSequences package exists, convert to appropriate format
  if(requireNamespace('arulesSequences', quietly = TRUE)){
    # build sequences in the required format
    library(arulesSequences)
    # prepare sequence objects
    s <- as(seqdata, 'transactions')
    # use cspade
    cs <- cspade(s, parameter = list(support = support, maxsize = maxsize))
    return(cs)
  } else {
    stop('Install arulesSequences for sequence mining (cspade)')
  }
}

# -----------------------------
# 6. Statistical significance helpers
# -----------------------------

# t-test (two groups) for numeric outcome
ttest_by_group <- function(df, outcome, group){
  fml <- as.formula(paste0(outcome, '~', group))
  res <- t.test(fml, data = df)
  broom::tidy(res)
}

# ANOVA for more than 2 groups
anova_by_group <- function(df, outcome, group){
  fml <- as.formula(paste0(outcome, '~', group))
  fit <- aov(fml, data = df)
  broom::tidy(fit)
}

# Chi-square test for two categorical variables
chi_sq <- function(df, var1, var2){
  tbl <- table(df[[var1]], df[[var2]])
  res <- chisq.test(tbl)
  list(tidy = broom::tidy(res), observed = res$observed, expected = res$expected)
}

# Permutation test (difference in means)
permutation_test_mean <- function(x, y, nperm = 5000, seed = 123){
  set.seed(seed)
  obs <- mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE)
  xy <- c(x, y)
  n <- length(x)
  perms <- replicate(nperm, {
    s <- sample(xy)
    mean(s[1:n], na.rm=TRUE) - mean(s[(n+1):length(xy)], na.rm = TRUE)
  })
  pval <- mean(abs(perms) >= abs(obs))
  list(observed = obs, pvalue = pval, perms = perms)
}

# -----------------------------
# 7. Clustering
# -----------------------------

cluster_kmeans <- function(df, numeric_cols = NULL, k = 3, scale_data = TRUE){
  if(is.null(numeric_cols)) numeric_cols <- names(df)[sapply(df, is.numeric)]
  mat <- df[, numeric_cols, drop = FALSE]
  if(scale_data) mat <- scale(mat)
  km <- kmeans(mat, centers = k)
  list(kmeans = km, cluster = km$cluster)
}

cluster_hclust <- function(df, numeric_cols = NULL, method = 'ward.D2', k = 3, scale_data = TRUE){
  if(is.null(numeric_cols)) numeric_cols <- names(df)[sapply(df, is.numeric)]
  mat <- df[, numeric_cols, drop = FALSE]
  if(scale_data) mat <- scale(mat)
  d <- dist(mat)
  h <- hclust(d, method = method)
  cutree(h, k = k)
}

cluster_dbscan <- function(df, numeric_cols = NULL, eps = 0.5, minPts = 5, scale_data = TRUE){
  if(is.null(numeric_cols)) numeric_cols <- names(df)[sapply(df, is.numeric)]
  mat <- as.matrix(df[, numeric_cols, drop = FALSE])
  if(scale_data) mat <- scale(mat)
  db <- dbscan::dbscan(mat, eps = eps, minPts = minPts)
  db$cluster
}

# -----------------------------
# 8. Utility: simple report generator
# -----------------------------

quick_report <- function(df, file = NULL){
  # Print diagnostics and run a few common analyses on the dataframe
  diag <- diagnose_data(df)
  cat('\nRunning common outlier detection on numeric cols...\n')
  nums <- names(df)[sapply(df, is.numeric)]
  outl <- lapply(df[nums], function(x) list(iqr = length(detect_outliers_iqr(x)), z = length(detect_outliers_zscore(x))))
  print(outl)
  cat('\nDetecting localization issues...\n')
  print(detect_localization_issues(df))
  if(!is.null(file)){
    saveRDS(list(diag = diag, outliers = outl), file = file)
    cat('Saved quick report to', file, '\n')
  }
  invisible(list(diag = diag, outliers = outl))
}
