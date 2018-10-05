FROM continuumio/miniconda3:latest

RUN conda install jupyterlab
RUN conda install -c conda-forge nodejs

RUN jupyter labextension install @jupyter-widgets/jupyterlab-manager

RUN conda install python=3.6

# See https://stackoverflow.com/questions/52582563/pip-install-attributeerror-distinfodistribution-dep-map
RUN conda install 'testpath<0.4'
RUN pip install --upgrade pip

RUN apt-get update && apt-get install -y build-essential

### Add this repo
ADD . /working/spikewidgets
WORKDIR /working/spikewidgets
RUN pip install .

### spikeinterface (jeremy branch)
RUN mkdir -p /working
WORKDIR /working
RUN git clone https://github.com/colehurwitz31/spikeinterface
WORKDIR /working/spikeinterface
RUN git checkout jeremy
RUN pip install .

WORKDIR /working/spikewidgets


