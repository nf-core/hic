FROM nfcore/base
LABEL authors="Nicolas Servant" \
      description="Docker image containing all requirements for nf-core/hic pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-hic-1.0dev/bin:$PATH
