#install.packages("devtools")
devtools::install_github("https://github.com/CostaLab/CrossTalkeR", build_vignettes = TRUE)                     
require(CrossTalkeR)  
vignette('CrossTalkeR')

# paths <- c('CTR' = system.file("extdata",
#                                "ctr_nils_bm_human_newformat.csv",
#                                package = "CrossTalkeR"),
#            'EXP' = system.file("extdata",
#                                "exp_nils_bm_human_newformat.csv",
#                                package = "CrossTalkeR"))

genes <- c("Crabp1", "Des", "Rgs5", "C1qb", "Pf4", "Eln", "Ogn", "Cd93", "Pecam1")

#output <- system.file("extdata", package = "CrossTalkeR")

# data <- generate_report(paths,
#                         genes,
#                         out_path=paste0(output,'/'),
#                         threshold=0,
#                         out_file = 'vignettes_example.html',
#                         output_fmt = "html_document",
#                         report = TRUE)


suppressPackageStartupMessages({require(CrossTalkeR)})
