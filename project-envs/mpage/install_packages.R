setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
script_dir <- dirname(script_path)

base_packages <- scan(file.path(script_dir, "CRAN_packages.txt"), what="", sep="\n")
bioc_packages <- scan(file.path(script_dir, "bioconductor_packages.txt"), what="", sep="\n")
github_packages <-  read.table(file.path(script_dir, "devtools_packages.txt"), 
                               header = FALSE, 
                               stringsAsFactors = FALSE)





tryCatch({
    install.packages(base_packages)
}, error = function(e) {
    message('CRAN R packages: ', conditionMessage(e))
    quit(status = 1)
})


tryCatch({
    BiocManager::install(bioc_packages)
}, error = function(e) {
    message('Bioconductor: ', conditionMessage(e))
    quit(status = 1)
})

tryCatch({
    lapply(github_packages$V1, function(pkg) devtools::install_github(pkg))
}, error = function(e) {
    message('github R packages: ', conditionMessage(e))
    quit(status = 1)
})

IRkernel::installspec(name = 'VSCODE_R', displayname = 'VSCODE_R', user = FALSE)