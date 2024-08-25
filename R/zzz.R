rmsk_profiler <- NULL
pybedtools <- NULL

.onLoad <- function(libname, pkgname) {
  path <- system.file("python", package = "rmskProfiler")
  rmsk_profiler <<- reticulate::import_from_path("rmsk_profiler", path = path, delay_load = TRUE)
  pybedtools <<- reticulate::import_from_path("pybedtools", path = path, delay_load = TRUE)
}
