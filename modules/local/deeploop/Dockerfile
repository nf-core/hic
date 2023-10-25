FROM python:3.5

RUN wget -q https://github.com/JinLabBioinfo/DeepLoop/archive/refs/heads/master.zip
RUN unzip master.zip

WORKDIR /DeepLoop-master
RUN pip3 install -r requirements.txt
