FROM registry.codeocean.com/codeocean/matlab:2017b-ubuntu16.04

ARG DEBIAN_FRONTEND=noninteractive
ARG MLM_LICENSE_FILE

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
      "build-essential=12.1ubuntu2" \
      "gcc=4:5.3.1-1ubuntu1" \
      "gpp=2.24-3" \
      "libnss3=2:3.28.4-0ubuntu0.16.04.6" \
      "libxcursor1=1:1.1.14-1ubuntu0.16.04.2" \
      "libxi6=2:1.7.6-1" \
      "libxtst6=2:1.2.2-1" \
    && rm -rf /var/lib/apt/lists/*
