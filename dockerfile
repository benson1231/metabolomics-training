FROM hidakakombu/metaboanalystr_4.0.0:1.2

RUN R -e 'BiocManager::install(c("ggrepel", "fitdistrplus", "RJSONIO"))'

RUN R -e 'devtools::install_github("xia-lab/OptiLCMS")'