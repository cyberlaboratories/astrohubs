FROM cyberhubs/corehub
MAINTAINER Falk Herwig

ENV DEBIAN_FRONTEND noninteractive

EXPOSE 8888

USER root

WORKDIR /root

ENV HOME /root

COPY apt_packages_mesahub.txt /tmp/apt_packages_mesahub.txt

RUN apt-get update && \
    apt-get -y install $(cat /tmp/apt_packages_mesahub.txt) && \
    apt-get autoremove --yes && apt-get clean all && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN pip install --upgrade pip && \
    pip2 install --upgrade pip && \
    pip3 install --upgrade pip 

RUN echo "Installing Python packages..."

ADD ./single_pip2_requirements_mesahub.txt /tmp/requirements/single_pip2_requirements_mesahub.txt
ADD ./single_pip3_requirements_mesahub.txt /tmp/requirements/single_pip3_requirements_mesahub.txt

RUN pip2 install -r /tmp/requirements/single_pip2_requirements_mesahub.txt
RUN pip3 install -r /tmp/requirements/single_pip3_requirements_mesahub.txt

RUN python2 -m pip install ipykernel && \
    python2 -m ipykernel install && python2 -m bash_kernel.install && \
    python3 -m pip install ipykernel && \
    python3 -m ipykernel install && python3 -m bash_kernel.install
    
COPY ./source/mesa.tar.gz  /tmp/mesa.tar.gz

RUN tar zxvf /tmp/mesa.tar.gz -C /tmp/

COPY jupyterhub-singleuser /usr/local/bin/jupyterhub-singleuser

COPY singleuser.sh /srv/singleuser/singleuser.sh

RUN  chmod -R 775 /home/user/notebooks

<<<<<<< HEAD
RUN (cd /tmp; wget http://www.astro.wisc.edu/~townsend/resource/download/mesasdk/mesasdk-x86_64-linux-20160129.tar.gz) &&\
    tar xvfz /tmp/mesa/mesasdk-x86_64-linux-20160129.tar.gz -C /opt &&\
=======
RUN tar xvfz /tmp/mesa/mesasdk-x86_64-linux-20160129.tar.gz -C /opt &&\
>>>>>>> f2556d8d8ebcafec88e48b511eb2d6d571efb16b
    chmod 775 /tmp/mesa/mesa_utils_install.sh &&\
    /tmp/mesa/mesa_utils_install.sh &&\
    rm -rf /tmp/mesa

RUN ([ ! -d /opt/hdf5-1.8.3 ] && (cd /tmp ; \
     wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.3/src/hdf5-1.8.3.tar.gz ;\
    tar -xzvf hdf5-1.8.3.tar.gz; cd hdf5-1.8.3 ; ./configure --prefix=/opt/hdf5-1.8.3;\
    make ; make install; cd ..; rm -rf hdf5*))

RUN ([ ! -d /opt/se-1.2 ] && (cd /tmp; \
    wget http://nugridstars.org/data-and-software/se-library-release/se-1.2.1.tar.gz ;\
    tar -xzvf se-1.2.1.tar.gz; cd se-1.2; \
    ./configure --prefix=/opt/se-1.2 --with-hdf5=/opt/hdf5-1.8.3;\
    make; make install; cd ..; rm -rf se*))

# Jupyter lab section

RUN jupyter serverextension enable --py jupyterlab --sys-prefix &&\
    jupyter labextension install @jupyterlab/hub-extension@0.6.0 

COPY .bash_aliases /home/user
COPY README_mesa.ipynb  /home/user/notebooks

RUN chown -R user:user /home/user && \
    chmod -R 775 /home/user/ && \
    chmod 775 /usr/local/bin/jupyterhub-singleuser && \
    chmod 775 /srv/singleuser/singleuser.sh

USER user
ENV HOME /home/user
WORKDIR /home/user

COPY ipython_config.py /etc/ipython/ipython_config.py

CMD ["/bin/bash", "/srv/singleuser/singleuser.sh"]

