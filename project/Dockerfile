# Usa Ubuntu 22.04 come base
FROM ubuntu:22.04

# Evita richieste interattive (es. timezone)
ENV DEBIAN_FRONTEND=noninteractive

# Installa R e librerie di sistema necessarie
RUN apt-get update && apt-get install -y \
    tzdata \
    r-base \
    build-essential \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libxt-dev \
    wget \
    git \
    && rm -rf /var/lib/apt/lists/*

# Installa BiocManager e pacchetti R
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')" \
    && R -e "BiocManager::install(c('GenomicRanges','IRanges','GenomeInfoDb','S4Vectors','rtracklayer','GEOquery','sva','DESeq2','apeglm'), ask=FALSE)" \
    && R -e "install.packages(c('ggplot2','dplyr','data.table','pheatmap'), repos='https://cloud.r-project.org')"

# Installa supporto per RMarkdown
RUN apt-get update && apt-get install -y pandoc \
    && R -e "install.packages(c('rmarkdown','knitr','tinytex'), repos='https://cloud.r-project.org')"

# Imposta la cartella di lavoro dentro il container
WORKDIR /home/project

# Comando di default
CMD [ "bash" ]
