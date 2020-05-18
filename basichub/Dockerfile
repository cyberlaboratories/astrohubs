FROM cyberhubs/corehub:1.2.4
MAINTAINER Falk Herwig

ENV DEBIAN_FRONTEND noninteractive

EXPOSE 8888

USER root

WORKDIR /root

ENV HOME /root

RUN echo "Installing Ubuntu apt-get packages..."

COPY apt_packages.txt   /tmp/apt.txt

RUN apt-get update && \
    apt-get -y install $(cat /tmp/apt.txt) && \
    apt-get autoremove --yes && apt-get clean all && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN pip install --upgrade pip && \
    pip3 install --upgrade pip 

RUN echo "Installing Python packages..."

ADD ./single_pip3_requirements_ppmstar.txt /tmp/requirements/single_pip3.txt

RUN pip3 install -r /tmp/requirements/single_pip3.txt 

RUN jupyter serverextension enable --py nbdime  && \
    jupyter nbextension install --py nbdime && \
    jupyter nbextension enable --py nbdime && \
    jupyter labextension install nbdime-jupyterlab

COPY jupyterhub-singleuser /usr/local/bin/jupyterhub-singleuser
COPY singleuser.sh /srv/singleuser/singleuser.sh

RUN chown -R user:user /home/user && \
    chmod -R 775 /home/user/ && \
    chmod 775 /usr/local/bin/jupyterhub-singleuser && \
    chmod 775 /srv/singleuser/singleuser.sh

USER user
ENV HOME /home/user
WORKDIR /home/user

CMD ["/bin/bash", "/srv/singleuser/singleuser.sh"]

RUN ln -s /user /home/user
RUN ln -s /data /home/user
