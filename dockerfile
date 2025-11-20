FROM hidakakombu/metaboanalystr_4.0.0:1.2

RUN R -e 'devtools::install_github("xia-lab/OptiLCMS")'

RUN R -e 'BiocManager::install(c("ggrepel", "fitdistrplus", "RJSONIO", "ellipse", "vegan", "factoextra", "pls", "som", "randomForest", "pheatmap"))'

RUN apt-get update && apt-get install -y curl