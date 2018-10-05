FROM magland/jp_proxy_widget:20180831

### Install conda packages
RUN conda install python=3.6

### Add this repo
ADD . /working
WORKDIR /working

RUN pip install .
