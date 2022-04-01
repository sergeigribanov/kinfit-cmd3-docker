FROM ssgribanov/hep-env:latest

MAINTAINER Sergei Gribanov <ssgribanov@gmail.com>

USER $USER
COPY yadisk.py /home/$USER/workdir/yadisk.py
COPY utils   /home/$USER/workdir/utils
COPY --chown=hep:hep notebooks /home/$USER/workdir/notebooks
ENV PYTHONPATH /home/$USER/workdir/utils:$PYTHONPATH
COPY install.sh install.sh
RUN sh install.sh
EXPOSE 8765

CMD [ "/bin/bash" ]
ENTRYPOINT jupyter notebook --no-browser  --port 8765