FROM jupyter/minimal-notebook AS base

# Install gammapy
RUN mamba install gcc jupyterlab "gammapy==1.2" ipykernel --yes
WORKDIR /gammapy-tools


# Clone V2DL3
WORKDIR /gammapy-tools
RUN git clone https://github.com/VERITAS-Observatory/V2DL3.git
WORKDIR /gammapy-tools/V2DL3

# Because its an env install of requirements
# Note the python version in the latest test throws issues with pytest
RUN grep -A 100 "dependencies:" environment-eventdisplay.yml | grep "-" | grep -v "python" | awk '{print $2}' | xargs mamba install --yes

# root_numpy throws issues too. Only VEGAS uses it
RUN mv setup.py _setup.py && grep -v "root_numpy" _setup.py > setup.py && pip install .


# Final product
FROM jupyter/minimal-notebook AS final

# Copy across the conda/mamba install
COPY --from=base /opt/conda /opt/conda

# Contents of repo
WORKDIR /gammapy-tools/tmp_build

# RUN gammapy download datasets
ENV GAMMAPY_DATA=/gammapy-tools/gammapy-datasets/1.2/
RUN mkdir -p $GAMMAPY_DATA
WORKDIR /gammapy-tools/
RUN gammapy download datasets

# Add package
ADD --chown=1000:100 . /gammapy-tools/tmp_build/gammapy-tools
WORKDIR /gammapy-tools/tmp_build/gammapy-tools



# RUN ls -lah
RUN pip install .
# RUN ./gammapy_tools/Hipparcos_MAG8_1997.dat $GAMMAPY_DATA/catalogs/
RUN cp /opt/conda/lib/python3.11/site-packages/gammapy_tools/Hipparcos_MAG8_1997.dat  $GAMMAPY_DATA/
RUN wget https://raw.githubusercontent.com/gammapy/gammapy/main/gammapy/datasets/map.py -O /opt/conda/lib/python3.11/site-packages/gammapy/datasets/map.py


USER root
RUN mkdir /local_data

# Clean up
RUN rm -r /gammapy-tools/tmp_build
USER jovyan
RUN mamba clean -a --yes
WORKDIR /local_data

# Keep alive for Docker ipykernel usage
RUN echo -e "#!/bin/bash\nwhile true; do sleep 5; done" >> /gammapy-tools/keep_alive.sh ; chmod a+x /gammapy-tools/keep_alive.sh 
CMD ["/gammapy-tools/keep_alive.sh"]