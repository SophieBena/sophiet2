FROM registry.gitlab.com/enki-portal/thermoengine:master
COPY start.sh start-notebook.sh start-singleuser.sh / usr/local/bin #buildkit
USER root
RUN chown -R ${NB_UID} ${HOME}
RUN pip install --no-cache-dir appmode
RUN jupyter nbextension enable --py --sys-prefix appmode
RUN jupyter serverextension enable --py --sys-prefix appmode
USER ${NB_USER}

#EOF
