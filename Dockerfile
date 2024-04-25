FROM ubuntu:22.04

RUN apt update && apt install -y build-essential clang wget
COPY . /sbva/

WORKDIR /sbva/
RUN wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz && tar xf eigen-3.4.0.tar.gz
RUN make
