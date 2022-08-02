FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:c57f-main

ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

RUN apt install -y software-properties-common
RUN apt remove r-base
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-key '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7'
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/debian bullseye-cran40/'

RUN apt update
RUN apt install r-recommended r-base-core r-base libcurl4-openssl-dev libssl-dev wget curl pandoc pandoc-citeproc --yes

RUN Rscript -e "install.packages('BiocManager')"
RUN Rscript -e "BiocManager::install('alevinQC')"
RUN Rscript -e "BiocManager::install('celda')"
RUN Rscript -e "BiocManager::install('zellkonverter')"

RUN curl -O \
  https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
  && mkdir /root/.conda \
  && bash Miniconda3-latest-Linux-x86_64.sh -b \
  && rm -f Miniconda3-latest-Linux-x86_64.sh

RUN conda install python=3.8
RUN conda install -c defaults \
  -c conda-forge \
  -c bioconda \
  -y -n base --debug \
  salmon alevin-fry==0.7.0 multiqc fastqc pyroe bedtools gffread


RUN python3 -m pip install --upgrade latch lgenome mygene requests openpyxl

COPY scDeepSort-v1.0-cpu.tar.gz /root/scDeepSort-v1.0-cpu.tar.gz

RUN echo 'deb http://deb.debian.org/debian testing main' > /etc/apt/sources.list.d/testing.list
RUN apt update
RUN apt install gcc-11 g++-11 -t testing -y
RUN rm /usr/lib/x86_64-linux-gnu/libstdc++.so.6
#RUN ls /usr/lib/x86_64-linux-gnu/libstdc++.so* && /bin/sleep 10
RUN ln /usr/lib/x86_64-linux-gnu/libstdc++.so.6.0.30 /usr/lib/x86_64-linux-gnu/libstdc++.so.6
RUN python3 -m pip install /root/scDeepSort-v1.0-cpu.tar.gz
ENV DGLBACKEND pytorch

RUN apt-get install libmagick++-dev -y
RUN Rscript -e "BiocManager::install('celda')"
RUN Rscript -e "package('celda)'"


COPY wf /root/wf
COPY scripts/qc.R /root/qc.R
COPY scripts/decontx.R /root/decontx.R
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
WORKDIR /root
