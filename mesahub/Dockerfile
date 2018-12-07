FROM cyberhubs/basichub

MAINTAINER Falk Herwig

ENV DEBIAN_FRONTEND noninteractive

EXPOSE 8888

USER root

WORKDIR /root

ENV HOME /root

COPY ./source/mesa.tar.gz  /tmp/mesa.tar.gz

RUN tar zxvf /tmp/mesa.tar.gz -C /tmp/

COPY singleuser.sh /srv/singleuser/singleuser.sh

RUN  chmod -R 775 /home/user/notebooks

RUN (cd /tmp/mesa; wget --user-agent=""  http://www.astro.wisc.edu/~townsend/resource/download/mesasdk/mesasdk-x86_64-linux-20160129.tar.gz) &&\
    tar xvfz /tmp/mesa/mesasdk-x86_64-linux-20160129.tar.gz -C /opt &&\
    chmod 775 /tmp/mesa/mesa_utils_install.sh &&\
    /tmp/mesa/mesa_utils_install.sh &&\
    rm -rf /tmp/mesa

RUN ([ ! -d /opt/hdf5-1.8.3 ] && (cd /tmp ; \
     wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.3/src/hdf5-1.8.3.tar.gz ;\
    tar -xzvf hdf5-1.8.3.tar.gz; cd hdf5-1.8.3 ; ./configure --prefix=/opt/hdf5-1.8.3;\
    make ; make install; cd ..; rm -rf hdf5*))

RUN (cd /tmp ;  wget https://www.open-mpi.org/software/ompi/v3.0/downloads/openmpi-3.0.0.tar.gz; \
    tar -xzvf openmpi-3.0.0.tar.gz ; cd openmpi-3.0.0;\
    mkdir /opt/openmpi-3.0.0; ./configure --prefix=/opt/openmpi-3.0.0;\
    make all install; cd ..; rm -rf openmpi*)

RUN echo Rerun from here: 2
RUN ([ ! -d /opt/se-1.2 ] && (cd /tmp; \
    git clone https://github.com/NuGrid/NuSE.git ;\
    cd NuSE/SE; \
    ./configure --prefix=/opt/se-1.2 --with-hdf5=/opt/hdf5-1.8.3;\
    make; make install; cd ../..; rm -rf NuSE))

COPY .bash_aliases /tmp
RUN cat /tmp/.bash_aliases >> /home/user/.bash_aliases

COPY README_mesa.ipynb  /home/user/notebooks

RUN chown -R user:user /home/user && \
    chmod -R 775 /home/user/ && \
    chmod 775 /usr/local/bin/jupyterhub-singleuser && \
    chmod 775 /srv/singleuser/singleuser.sh

USER user
ENV HOME /home/user
WORKDIR /home/user

CMD ["/bin/bash", "/srv/singleuser/singleuser.sh"]

RUN ln -s /home/user/notebooks/README_mesa.ipynb /home/user
