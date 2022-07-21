FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:9a7d-main

ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"


RUN apt install -y software-properties-common
RUN apt remove r-base
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-key '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7'
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/debian buster-cran40/'

RUN apt update
RUN apt install r-recommended r-base-core r-base libcurl4-openssl-dev libssl-dev wget curl pandoc pandoc-citeproc --yes

RUN Rscript -e "install.packages('BiocManager')"
RUN Rscript -e "BiocManager::install('alevinQC')"
RUN Rscript -e "library(alevinQC)"

RUN curl -O \
  https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
  && mkdir /root/.conda \
  && bash Miniconda3-latest-Linux-x86_64.sh -b \
  && rm -f Miniconda3-latest-Linux-x86_64.sh

RUN conda install -c defaults \
  -c conda-forge \
  -c bioconda \
  -y -n base --debug \
  salmon alevin-fry multiqc fastqc pyroe bedtools gffread

RUN python3 -m pip install latch==1.9.0 lgenome mygene
RUN pip install --upgrade requests
COPY wf /root/wf
COPY scripts/qc.R /root/qc.R
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
WORKDIR /root
