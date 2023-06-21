# FROM rocker/shiny-verse:latest
FROM rocker/shiny-verse:4.0.0

RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libglpk40 \
    libudunits2-dev \
    libproj-dev \
    libgdal-dev \
    libgeos-dev \
    libnode-dev \
    libzip-dev \
    bzip2

## Install R libraries
RUN R -e "install.packages('shiny', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('shinydashboard', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('remotes', repos='http://cran.rstudio.com/')"

RUN R -e "install.packages(c('ade4',  'adegenet',  'ape',  'BiocManager',  'castor', \
          'colourpicker',  'data.table',  'dplyr',  'DT', 'geosphere', \
          'ggplot2', 'ggtree',  'ggtree',  'globe4r',  'hash',  'htmltools', \
          'htmlwidgets',  'igraph', 'import',  'knitr',  'leaflet', \
          'magrittr',  'markdown',  'network', 'phangorn', 'plotly',  'plyr',  \
          'randomcoloR',  'rbokeh',  'readr',  'rhandsontable',  \
          'rmarkdown',  'shinycssloaders',  'shinyjqui',  \
          'shinythemes',  'shinyWidgets',  'stringr',  'tibble',  \
          'visNetwork',  'webshot'))"

RUN R -e "install.packages('bslib', repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('treeio')"

## Install StrainHub
RUN R -e "remotes::install_github('colbyford/strainhub', subdir='pkg', dependencies=TRUE)"
# RUN R -e "remotes::install_github('colbyford/strainhub', subdir='pkg', dependencies=FALSE)"


## Copy StrainHub application files
COPY /app /srv/shiny-server/
COPY /data /srv/shiny-server/data


## Open the port for Shiny
EXPOSE 3838

## Change permissions
RUN sudo chown -R shiny:shiny /srv/shiny-server

## Fix PNG export functionality with webshot
RUN R -e "webshot::install_phantomjs()"

## Start Shiny Server
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/', host = '0.0.0.0', port = 3838)"]