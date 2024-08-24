#' Install the Python dependencies for rmskProfiler
#'
#' This function simply wraps reticulate::py_install() for installing the needed
#' python packages into the default environment named "r-rmskProfiler". Depending
#' on your system setup it is likely that you'll need to pass additional
#' arguments. See the example how to install using a mambaforge env.
#'
#' @param envname The name, or full path, of the environment in which Python
#' packages are to be installed. The default is "r-rmskProfiler", probably don't
#' change that unless you have a good reason. When NULL, the active environment
#' as set by the RETICULATE_PYTHON_ENV variable will be used; if that is unset,
#' then the r-reticulate environment will be used.
#' @param method Installation method. By default, "auto" automatically finds a
#' method that will work in the local environment. Change the default to force
#' a specific installation method. Note that the "virtualenv" method is not
#' available on Windows.
#' @param ... Additional arguments passed to reticulate::py_install()
#' @export
#' @examples
#' \dontrun{
#'
#' # Install requirements using mambaforge and specifying bioconda channel
#' install_rmskProfiler(
#'   method = "conda",
#'   conda = "/home/gennaro/mambaforge/condabin/conda",
#'   channel = "bioconda"
#'   )
#' }
install_rmskProfiler <- function(envname = "r-rmskProfiler", method = "auto", ...) {
  reticulate::py_install(
    packages = c("pybedtools"),
    python_version = "3.10",
    envname = envname,
    method = method,
    ...
    )
}
