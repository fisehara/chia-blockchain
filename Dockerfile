FROM python:3.8-slim-buster
ARG BUILD_BRANCH=master

WORKDIR /usr/src

RUN apt-get update
RUN apt-get install --assume-yes build-essential git cmake libssl-dev libffi-dev libgmp3-dev python3-virtualenv --no-install-recommends

RUN git clone https://github.com/Chia-Network/chia-blockchain.git .
RUN git checkout $BUILD_BRANCH

RUN python3 -m venv .venv
RUN . .venv/bin/activate && pip install wheel
RUN . .venv/bin/activate && pip install -e .
RUN . .venv/bin/activate && pip install -r requirements.txt
RUN . .venv/bin/activate && python -m scripts.regenerate_keys