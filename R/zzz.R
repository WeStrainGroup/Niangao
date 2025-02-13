#' @title Package Initialization
#'
#' @description
#' This function is called when the package is loaded. It checks and installs
#' the required dependencies (`sangerseqR` and `Biostrings`) if they are not
#' already installed. It is automatically executed when the package is loaded.
#'
#' @param libname A character string representing the library name. This argument
#'   is passed to `.onLoad` by the R environment when the package is loaded.
#' @param pkgname A character string representing the package name. This argument
#'   is passed to `.onLoad` by the R environment when the package is loaded.
#'
#' @details
#' The function first checks if `BiocManager` is installed. If it is not, it will
#' install it. Then, it checks for the presence of the required dependencies,
#' namely `sangerseqR` and `Biostrings`. If any of these packages are missing,
#' it installs them using `BiocManager`.
#'
#' @seealso
#' \code{\link{check_and_install_dependencies}} for the package installation logic.
#'
#' @import BiocManager
.onLoad <- function(libname, pkgname) {
  check_and_install_dependencies()
}

#' @title Check and Install Missing Dependencies
#'
#' @description
#' This function checks whether the required dependencies (`sangerseqR` and `Biostrings`)
#' are installed on the user's system. If any of the packages are missing, it installs them
#' using `BiocManager`.
#'
#' @details
#' The function checks for the presence of the `BiocManager` package, and if it is not
#' installed, it installs it first. Then, it checks for the presence of `sangerseqR` and
#' `Biostrings`. If either of these packages is missing, the function installs them via
#' `BiocManager::install()`.
#'
#' @import BiocManager
#' @importFrom utils install.packages
check_and_install_dependencies <- function() {
  required_packages <- c("sangerseqR", "Biostrings")

  # Check and install BiocManager if not installed
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }

  # Check and install required packages if not installed
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Installing missing package:", pkg))
      BiocManager::install(pkg)
    }
  }
}
