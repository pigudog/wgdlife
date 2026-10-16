# R 4.4.3; the only package the convergence test needs (data.table); `parallel` ships with R.
if (!requireNamespace("data.table", quietly = TRUE) || packageVersion("data.table") != "1.18.6.1") {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = "https://cloud.r-project.org")
  remotes::install_version("data.table", version = "1.18.6.1", repos = "https://cloud.r-project.org")
}
cat("data.table", as.character(packageVersion("data.table")), "\n")
