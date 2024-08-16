FROM openanalytics/r-shiny

# Install Additional R Packages (if needed)
RUN R -e "install.packages(c('devtools', 'BiocManager', 'DT', 'tidyverse', 'janitor', 'ggsci', 'ggthemes', 'here','Matrix', 'readxl','RColorBrewer', 'ggheatmap','vegan', 'corto', 'BH', 'plotly','patchwork', 'flexdashboard', 'heatmaply'), dependencies = TRUE, repos='https://cloud.r-project.org')"

# Install Bioconductor Packages (if needed)
RUN R -e "BiocManager::install(c('clusterProfiler', 'org.Hs.eg.db', 'cytolib'))"

# Copy App Directory (replace 'VaxGO_Tool' with your actual directory name)
COPY . /srv/shiny-server

# Expose Shiny Server Port
EXPOSE 3838

# Start Shiny Server with your App
CMD ["R", "-e", "rmarkdown::run('VaxGO_Tool.Rmd', shiny_args = list(host = '0.0.0.0', port = 3838))"]

