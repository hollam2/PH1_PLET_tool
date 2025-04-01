
# List of required packages
required_packages <- c("reticulate", "testthat")

# Install missing packages
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
  }
}

# Install necessary packages
install_if_missing(required_packages)