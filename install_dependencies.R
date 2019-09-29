# PEXO dependencies
list.of.packages <- c("gmp", "orthopolynom", "pracma", "cubature", "HyperbolicDist", "magicaxis", "numDeriv", "optparse")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if (length(new.packages)) {
    print("Installing packages:")
    print(new.packages)
    install.packages(new.packages, repos = "http://cran.us.r-project.org")
} else {
    print("All packages are installed.")
}
