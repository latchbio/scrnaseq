FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:9a7d-main

ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

RUN apt-get update && apt-get install --yes \
  libcurl4-openssl-dev \
  libxml2-dev \
  libssl-dev \
  curl \
  software-properties-common \
  dirmngr

RUN curl -O \
  https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
  && mkdir /root/.conda \
  && bash Miniconda3-latest-Linux-x86_64.sh -b \
  && rm -f Miniconda3-latest-Linux-x86_64.sh

RUN conda install -c defaults \
  -c conda-forge \
  -c bioconda \
  -y -n base --debug \
  salmon alevin-fry multiqc fastqc pyroe bedtools

# STOP HERE:
# The following lines are needed to ensure your build environement works
# correctly with latch.
RUN python3 -m pip install --upgrade latch lgenome requests
COPY wf /root/wf
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
WORKDIR /root
CMD ["/bin/sleep", "1000"]
