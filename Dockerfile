FROM nfcore/base:1.9
LABEL authors="Nicolas Servant" \
      description="Docker image containing all software requirements for the nf-core/hic pipeline"

## Install gcc for pip iced install
RUN apt-get update && apt-get install -y gcc g++ && apt-get clean -y

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-hic-1.2.0/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-hic-1.2.0 > nf-core-hic-1.2.0.yml

