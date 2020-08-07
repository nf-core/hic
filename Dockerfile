FROM nfcore/base:1.10.2

LABEL authors="Nicolas Servant" \
      description="Docker image containing all software requirements for the nf-core/hic pipeline"

## Install gcc for pip iced install
RUN apt-get update && apt-get install -y gcc g++ && apt-get clean -y

COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-hic-1.2.1/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-hic-1.2.1 > nf-core-hic-1.2.1.yml

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron
