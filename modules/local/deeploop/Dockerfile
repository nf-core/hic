FROM python:3.6

RUN apt update && apt install -y liblist-moreutils-perl r-base r-base-dev

RUN wget -q https://github.com/JinLabBioinfo/DeepLoop/archive/refs/heads/master.zip
RUN unzip master.zip

WORKDIR /DeepLoop-master
RUN pip3 install -r requirements.txt

ENV PATH /DeepLoop-master/prediction:$PATH

RUN wget --no-check-certificate https://hiview.case.edu/ssz20/tmp.HiCorr.ref/DeepLoop_models.tar.gz
RUN tar -xvf DeepLoop_models.tar.gz

RUN chmod +x lib/*.pl lib/*.r
