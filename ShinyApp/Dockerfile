# Use rocker/shiny as the base image
FROM rocker/shiny:latest

# Set the working directory inside the container
WORKDIR /app

# Copy the application files into the container
COPY . /app

# Install system dependencies (if any)
# Example:
# RUN apt-get update && apt-get install -y --no-install-recommends libxml2-dev

# Install R package dependencies
RUN R -e "install.packages(c('shiny', 'Matrix', 'circlize', 'RColorBrewer', 'org.Hs.eg.db', 'plotly', 'clusterProfiler', 'janitor', 'ggsci', 'patchwork', 'ggthemes', 'here', 'readxl', 'scales', 'tidyverse', 'vegan', 'corto', 'DT', 'BH', 'cytolib', 'heatmaply', 'BiocManager', 'devtools'), repos='https://cran.r-project.org')"
RUN R -e "if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install(c('GEOquery', 'DESeq2', 'biomaRt','celldex','GSVA', 'sva','AnnotationDbi', 'ComplexHeatmap', 'InteractiveComplexHeatmap', 'ImmuneSpaceR'), ask = FALSE, update = FALSE)"

# Expose port 3838 for Shiny
EXPOSE 3838

# Command to run the Shiny application
CMD ["R", "-e", "shiny::runApp('/app', host = '0.0.0.0', port = 3838)"]