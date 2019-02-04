FROM cyberhubs/basichub

MAINTAINER Falk Herwig

ENV DEBIAN_FRONTEND noninteractive

EXPOSE 8888

USER root
WORKDIR /root
ENV HOME /root


COPY apt_packages_mesahub.txt /tmp/apt_packages.txt

RUN apt-get update && apt-get --no-install-recommends install --yes $(cat /tmp/apt_packages.txt)
RUN apt-get autoremove --yes && apt-get clean all
RUN rm -rf /var/lib/apt/lists/*  /var/tmp/*

RUN echo Installed apt packages ..........

RUN ([ ! -d /opt/hdf5-1.8.3 ] && (cd /tmp ; \
     wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.3/src/hdf5-1.8.3.tar.gz ;\
    tar -xzvf hdf5-1.8.3.tar.gz; cd hdf5-1.8.3 ; ./configure --prefix=/opt/hdf5-1.8.3;\
    make ; make install; cd ..; rm -rf hdf5*))

RUN echo Installed hdf5 ..........

RUN (cd /tmp ;  wget https://www.open-mpi.org/software/ompi/v3.0/downloads/openmpi-3.0.0.tar.gz; \
    tar -xzvf openmpi-3.0.0.tar.gz ; cd openmpi-3.0.0;\
    mkdir /opt/openmpi-3.0.0; ./configure --prefix=/opt/openmpi-3.0.0;\
    make all install; cd ..; rm -rf openmpi*)

RUN echo Installed openmpi ..........

RUN ([ ! -d /opt/se-1.2 ] && (cd /tmp; \
    git clone https://github.com/NuGrid/NuSE.git ;\
    cd NuSE/SE; \
    ./configure --prefix=/opt/se-1.2 --with-hdf5=/opt/hdf5-1.8.3;\
    make; make install; cd ../..; rm -rf NuSE))

RUN echo Installed nuse ..........

COPY singleuser.sh /srv/singleuser/singleuser.sh
RUN  chmod -R 775 /home/user/notebooks

COPY .bash_aliases /tmp
RUN cat /tmp/.bash_aliases >> /home/user/.bash_aliases

RUN chown -R user:user /home/user && \
    chmod -R 775 /home/user/ && \
    chmod 775 /usr/local/bin/jupyterhub-singleuser && \
    chmod 775 /srv/singleuser/singleuser.sh

COPY ./source/mesa.tgz  /tmp/mesa.tgz
RUN tar zxvf /tmp/mesa.tgz -C /tmp/ &&\
    chmod 775 /tmp/mesa/mesa_utils_install.sh  &&\
    /tmp/mesa/mesa_utils_install.sh &&\
    rm -rf /tmp/mesa*

RUN git clone https://github.com/NuGrid/NuPyCEE.git /opt/NuPyCEE && \
    git clone https://github.com/becot85/JINAPyCEE.git /opt/JINAPyCEE && \
    cd /opt/JINAPyCEE && python2 setup.py develop  && python3 setup.py develop
ENV PYTHONPATH /usr/local/bin/:/opt/NuPyCEE/:/opt

USER user:user
ENV HOME /home/user
WORKDIR /home/user

COPY mesasdk-x86_64-linux-20180822.tar.gz /home/user/
RUN tar xvfz /home/user/mesasdk-x86_64-linux-20180822.tar.gz -C /home/user
RUN rm -f /home/user/mesasdk-x86_64-linux-20180822.tar.gz

RUN (cd /home/user/notebooks;git clone https://github.com/NuGrid/wendi-examples.git;\
     ln -s /home/user/notebooks/wendi-examples /home/user)

CMD ["/bin/bash", "/srv/singleuser/singleuser.sh"]

COPY README_mesa.ipynb  /home/user/notebooks
COPY README_wendi3.ipynb  /home/user/notebooks
RUN ln -s /home/user/notebooks/README_mesa.ipynb /home/user &&\
    ln -s /home/user/notebooks/README_wendi3.ipynb /home/user
