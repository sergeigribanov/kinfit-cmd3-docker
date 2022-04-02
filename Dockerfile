FROM ssgribanov/hep-env:latest

MAINTAINER Sergei Gribanov <ssgribanov@gmail.com>

USER $USER
COPY utils   /home/$USER/workdir/utils
COPY --chown=$USER:$USER notebooks /home/$USER/workdir/notebooks
COPY --chown=$USER:$USER download_data.ipynb /home/$USER/workdir/download_data.ipynb
ENV PYTHONPATH /home/$USER/workdir/utils:$PYTHONPATH
COPY install.sh install.sh
RUN sh install.sh
USER 0
RUN mkdir /var/kinfit && chown $USER:$USER /var/kinfit
USER $USER
EXPOSE 8765

CMD [ "/bin/bash" ]
ENTRYPOINT jupyter notebook --no-browser  --port 8765