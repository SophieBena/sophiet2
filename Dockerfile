FROM registry.gitlab.com/enki-portal/thermoengine:master
COPY . ${HOME}
USER root
RUN chown -R ${NB_UID} ${HOME}
RUN pip install --no-cache-dir appmode
RUN jupyter nbextension enable --py --sys-prefix appmode
RUN jupyter serverextension enable --py --sys-prefix appmode
USER ${NB_USER}

FROM jupyter/base-notebook

# copy the notebook file to the container
COPY Ol-Opx-SplV14.ipynb /app/

# set the working directory
WORKDIR /app

# start the notebook
CMD ["jupyter", "notebook", "Ol-Opx-SplV14.ipynb"]

#EOF
