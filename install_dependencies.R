# PEXO dependencies
#list.of.packages <- c("orthopolynom", "pracma", "optparse")
list.of.packages <- c('optparse','orthopolynom','pracma','foreach','compiler','foreach','doMC','MASS')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if (length(new.packages)) {
    print("Installing packages:")
    print(new.packages)
    install.packages(new.packages, dependencies=TRUE, repos='http://cran.rstudio.com/')
} else {
    print("All packages are installed.")
}
