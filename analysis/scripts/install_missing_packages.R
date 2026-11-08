# install_missing_packages.R
# Purpose: install Mr.robust and MVMR (missing from baseline package audit) from
#          Posit Package Manager (P3M) binary repo for R 4.4.x on Windows.
# Usage:   Rscript analysis/scripts/install_missing_packages.R
# Note:    fallback sources (CRAN Archive / GitHub) handled manually if this fails.
options(repos = c(P3M = "https://packagemanager.posit.co/cran/latest"))
install.packages(c("Mr.robust", "MVMR"))
cat("--- verification ---\n")
for (p in c("Mr.robust", "MVMR")) {
  cat(p, "installed:", requireNamespace(p, quietly = TRUE), "\n")
}
cat("--- session ---\n")
cat(R.version.string, "\n")
