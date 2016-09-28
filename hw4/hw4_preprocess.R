########################################
# Code from James
########################################
# read svmlight format matrix for classification problem
# each row is in the following format
# label feature1:val1 feature2:val2 ... featureK:valK
# assumes one-indexing for the column indices
# This function returns a list with the labels as a vector
# and the features stored as a sparse Matrix (or a simple triplet matrix)
# It requires the Matrix and readr packages.
read_svmlight_class = function(myfile, format='sparseMatrix', num_cols = NULL) {
  require(Matrix)
  require(readr)

  raw_x = read_lines(myfile)
  x = strsplit(raw_x, ' ', fixed=TRUE)
  x = lapply(x, function(y) strsplit(y, ':', fixed=TRUE))
  l = lapply(x, function(y) as.numeric(unlist(y)))
  label = as.integer(lapply(l, function(x) x[1]))
  num_rows = length(label)
  features = lapply(l, function(x) tail(x,-1L))
  row_length = as.integer(lapply(features, function(x) length(x)/2))
  features = unlist(features)
  i = rep.int(seq_len(num_rows), row_length)
  j = features[seq.int(1, length(features), by = 2)] + 1
  v = features[seq.int(2, length(features), by = 2)]

  if(missing(num_cols)) {
    num_cols = max(j)
  }
  m = sparseMatrix(i=i, j=j, x=v, dims=c(num_rows, num_cols))

  list(labels=label, features=m)
}

# Where are the files stored?
base_dir = "~/Desktop/sds385/url_svmlight/"
svm_files = dir(base_dir, pattern = "*.svm")

# Loop through the files and create a list of objects
X_list = list()
y_list = list()
for(i in seq_along(svm_files)) {
  myfile = svm_files[i]
  cat(paste0("Reading file ", i, ": ", myfile, "\n"))
  D = read_svmlight_class(paste0(base_dir, myfile), num_cols = 3231962)
  X_list[[i]] = D$features
  y_list[[i]] = D$labels
}

# Assemble one matrix of features/vector of responses (do.call very handy here, although not super efficient
X = do.call(rBind, X_list)  # rBind, not rbind, for sparse matrices
y = do.call(c, y_list)
y = 0 + {y==1}

# Save as serialized (binary) files for much faster read-in next time
saveRDS(X, file='url_X.rds')
saveRDS(y, file='url_y.rds')