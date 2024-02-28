FROM mercury/azimuth_dsb:latest

RUN apt install cmake -y
RUN R --no-echo -e "install.packages('nloptr', dependencies = TRUE)"
RUN R --no-echo -e "install.packages('ggpubr')"
RUN R --no-echo -e "install.packages('viridis')"

