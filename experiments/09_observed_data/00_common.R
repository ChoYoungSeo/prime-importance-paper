load_observed_data <- function(id, cache_dir = "data/real_data_cache") {
  clean <- function(X, y, label) {
    X <- as.matrix(X)
    storage.mode(X) <- "double"
    y <- as.numeric(y)
    keep <- is.finite(y) & apply(X, 1L, function(row) all(is.finite(row)))
    X <- X[keep, , drop = FALSE]
    y <- y[keep]
    nonconstant <- apply(X, 2L, stats::sd) > 0
    list(
      id = id, label = label,
      X = X[, nonconstant, drop = FALSE], y = y
    )
  }

  if (id == "concrete") {
    data <- utils::read.csv(file.path(cache_dir, "concrete.csv"))
    return(clean(
      data[, setdiff(names(data), "strength")], data$strength,
      "Concrete Compressive Strength"
    ))
  }
  if (id == "airfoil") {
    data <- utils::read.table(
      file.path(cache_dir, "airfoil_self_noise.dat"), header = FALSE
    )
    names(data) <- c(
      "Frequency", "Angle", "Chord_Length", "Velocity",
      "Displacement", "Sound_Pressure"
    )
    return(clean(
      data[, 1:5], data$Sound_Pressure, "Airfoil Self-Noise"
    ))
  }
  if (id == "auto_mpg") {
    data <- utils::read.table(
      file.path(cache_dir, "auto-mpg.data"),
      header = FALSE, na.strings = "?", stringsAsFactors = FALSE,
      quote = "\""
    )
    names(data) <- c(
      "mpg", "cylinders", "displacement", "horsepower", "weight",
      "acceleration", "model_year", "origin", "car_name"
    )
    data <- stats::na.omit(data[, 1:8])
    data$origin_2 <- as.numeric(data$origin == 2)
    data$origin_3 <- as.numeric(data$origin == 3)
    data$origin <- NULL
    return(clean(
      data[, setdiff(names(data), "mpg")], data$mpg, "Auto MPG"
    ))
  }
  stop("Unknown observed dataset: ", id, call. = FALSE)
}

inject_permuted_controls <- function(X, number, seed) {
  X <- as.matrix(X)
  set.seed(seed)
  n <- nrow(X)
  parents <- sample.int(ncol(X), number, replace = TRUE)
  controls <- vapply(
    parents,
    function(parent) X[sample.int(n), parent],
    numeric(n)
  )
  colnames(controls) <- paste0("Permuted_", seq_len(number))
  list(
    X = cbind(X, controls),
    original = seq_len(ncol(X)),
    controls = ncol(X) + seq_len(number)
  )
}

atomic_save_observed <- function(object, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(object, temporary)
  if (!file.rename(temporary, path)) {
    unlink(temporary)
    stop("Could not save checkpoint: ", path, call. = FALSE)
  }
}
