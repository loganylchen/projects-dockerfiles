message('install something')
if (!requireNamespace('devtools', quietly = TRUE))
    install.packages('devtools')

devtools::install_github('kevinblighe/PCAtools')

