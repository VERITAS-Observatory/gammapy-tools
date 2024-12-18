FROM python:3.12-slim-bookworm AS base


RUN apt update && \
    apt upgrade -y && \
    apt install git curl bash wget gcc -y

RUN python3 -m venv /opt/gammapy_tools && \
    . /opt/gammapy_tools/bin/activate && \
    pip install --no-cache iminuit cmasher pip papermill matplotlib \
        uproot "jupyterlab==4.0.12"  notebook ipykernel ipython ipywidgets "gammapy==1.3"
    
RUN mkdir /gamma-tools && cd /gamma-tools && \
    git clone  --depth 1 --branch v0.6.0 https://github.com/VERITAS-Observatory/V2DL3.git && \
    cd V2DL3 && \
    grep -A 100 "dependencies:" environment-eventdisplay.yml | grep "-" | grep -v "python" | awk '{print $2}' | xargs pip install && \
    mv setup.py _setup.py && grep -v "root_numpy" _setup.py > setup.py && pip install . && \
    cd /gamma-tools/ && rm -r /gamma-tools/V2DL3

WORKDIR /gamma-tools/
RUN /opt/gammapy_tools/bin/gammapy download datasets

COPY . /gamma-tools/tmp_build
RUN cd /gamma-tools/tmp_build ; pip install . ;  cd /gamma-tools ; rm -r /gamma-tools/tmp_build && \
    python -m pip cache purge


# FROM python:3.12-slim-bookworm AS final

# COPY --from=base /opt/gammapy_tools /opt/gammapy_tools
# COPY --from=base /gamma-tools /gamma-tools

# ENV JPORT=8000
# ENV GAMMAPY_DATA=/gamma-tools/gammapy-datasets/1.3/