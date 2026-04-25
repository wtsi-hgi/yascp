FROM mercury/azimuth_dsb

RUN apt install cmake -y
RUN R --no-echo -e "install.packages('nloptr', dependencies = TRUE)"
RUN R --no-echo -e "install.packages('ggpubr')"
RUN R --no-echo -e "install.packages('viridis')"
RUN R --no-echo -e "remotes::install_version('Matrix', '1.6.1');"
RUN R --no-echo -e "install.packages('Rools');"
RUN R --no-echo -e "install.packages('irlba',type='source', dependencies=T)"
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
RUN apt-get update
RUN apt-get install r-base r-base-dev -y
RUN R --version
RUN apt-get update && apt-get install -y libglpk-dev
RUN apt install gfortran
RUN R --no-echo -e 'install.packages("igraph")'
RUN R --no-echo -e "install.packages('tidygraph')"
RUN R --no-echo -e "devtools::install_github('immunomind/immunarch', ref='dev')"
