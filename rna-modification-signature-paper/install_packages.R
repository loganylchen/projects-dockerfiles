setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
script_dir <- dirname(script_path)


github_packages <-  read.table(file.path(script_dir, "devtools_packages.txt"), 
                               header = FALSE, 
                               stringsAsFactors = FALSE)



# BiocManager::install("pd.hg.u133.plus.2")
BiocManager::install("oligo")

tryCatch({
    lapply(github_packages$V1, function(pkg) devtools::install_github(pkg))
}, error = function(e) {
    message('github R packages: ', conditionMessage(e))
    quit(status = 1)
})



IRkernel::installspec(name = 'VSCODE_R', displayname = 'VSCODE_R', user = FALSE)
