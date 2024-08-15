FROM rocker/r-ver:latest

LABEL maintainer="wasimvacinas@gmail.com"

# Instalar pacotes do sistema
RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    && rm -rf /var/lib/apt/lists/*

# Instalar Shiny
RUN R -e "install.packages('shiny', repos='https://cloud.r-project.org')"

# Instale pacotes R de CRAN
RUN R -e "install.packages(c('devtools', 'BiocManager', 'DT', 'janitor', 'ggsci', 'ggthemes', 'here', 'readxl', 'scales', 'vegan', 'corto', 'plotly', 'flexdashboard', 'tidyverse', 'heatmaply'))"

# Instale pacotes do Bioconductor
RUN R -e "BiocManager::install(c('DESeq2', 'GEOquery', 'biomaRt', 'clusterProfiler', 'ComplexHeatmap', 'InteractiveComplexHeatmap', 'org.Hs.eg.db', 'GSVA', 'sva', 'EnhancedVolcano', 'cytolib', 'ImmuneSpaceR'))"

# Copie o conteúdo do diretório atual para o container
COPY . /srv/shiny-server/

# Exponha a porta padrão do Shiny
EXPOSE 3838

# Defina o ponto de entrada para o Shiny
CMD ["R", "-e", "rmarkdown::run('/srv/shiny-server/VaxGO_Tool.Rmd', shiny_args = list(host = '0.0.0.0', port = 3838))"]
