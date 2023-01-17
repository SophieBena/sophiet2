FROM registry.gitlab.com/enki-portal/thermoengine:master
COPY inital-condarc / opt/conda/.condarc # buildkit
USER root
RUN chown -R ${NB_UID} ${HOME}
RUN pip install --no-cache-dir appmode
RUN jupyter nbextension enable --py --sys-prefix appmode
RUN jupyter serverextension enable --py --sys-prefix appmode
USER ${NB_USER}

#EOF
