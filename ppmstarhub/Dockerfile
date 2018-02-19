FROM cyberhubs/corehub
MAINTAINER Falk Herwig

ENV DEBIAN_FRONTEND noninteractive

EXPOSE 8888

USER root

WORKDIR /root

ENV HOME /root

COPY apt_packages_ppmstar.txt   /tmp/apt.txt

RUN apt-get update && \
    apt-get -y install $(cat /tmp/apt.txt) && \
    apt-get autoremove --yes && apt-get clean all && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN pip install --upgrade pip && \
    pip2 install --upgrade pip && \
    pip3 install --upgrade pip 

RUN echo "Installing Python packages..."

ADD ./single_pip2_requirements_ppmstar.txt /tmp/requirements/single_pip2.txt
ADD ./single_pip3_requirements_ppmstar.txt /tmp/requirements/single_pip3.txt



RUN pip2 install -r /tmp/requirements/single_pip2.txt  
RUN pip3 install -r /tmp/requirements/single_pip3.txt 

RUN python2 -m pip install ipykernel && \
    python2 -m ipykernel install && python2 -m bash_kernel.install && \
    python3 -m pip install ipykernel && \
    python3 -m ipykernel install && python3 -m bash_kernel.install
    

COPY jupyterhub-singleuser /usr/local/bin/jupyterhub-singleuser
COPY singleuser.sh /srv/singleuser/singleuser.sh
RUN  chmod -R 775 /home/user/notebooks

# moments and e3d stuff

COPY moments.tar.gz /tmp/moments.tar.gz
COPY src_read_ppm /tmp/src_read_ppm
COPY src_e3d /tmp/src_e3d

RUN make -C /tmp/src_read_ppm gfort && \
    mv /tmp/src_read_ppm/read_ppm /bin/read_ppm && \
    make -C /tmp/src_e3d gfort && \
    mv /tmp/src_e3d/e3d /bin/e3d
    
RUN pip2 install --upgrade pip && \
    pip2 install /tmp/moments.tar.gz && \
    pip3 install --upgrade pip && \
    pip3 install /tmp/moments.tar.gz && \
    rm -r /tmp/moments.tar.gz 

# Jupyter lab section

RUN jupyter serverextension enable --py jupyterlab --sys-prefix &&\
    jupyter labextension install @jupyterlab/hub-extension@0.6.0 

COPY .bash_aliases /home/user
COPY README_ppmstar.ipynb  /home/user/notebooks

RUN mkdir /home/user/notebooks/PPMstar && \
    git clone https://github.com/PPMstar/PPMnotebooks.git /home/user/notebooks/PPMstar/PPMnotebooks

#ADD PPMnotebooks /home/user/notebooks/PPMstar

RUN chown -R user:user /home/user && \
    chmod -R 775 /home/user/ && \
    chmod 775 /usr/local/bin/jupyterhub-singleuser && \
    chmod 775 /srv/singleuser/singleuser.sh

USER user
ENV HOME /home/user
WORKDIR /home/user

COPY ipython_config.py /etc/ipython/ipython_config.py

CMD ["/bin/bash", "/srv/singleuser/singleuser.sh"]

RUN ln -s /home/user/notebooks/README_ppmstar.ipynb /home/user

RUN ln -s /user/dataA_rpod4 /home/user/notebooks
RUN ln -s /user/dataA_rpod4 /home/user
RUN ln -s /data/cedar_projects /home/user/notebooks
RUN ln -s /data/cedar_projects /home/user
RUN ln -s /data/ppm_rpod2 /home/user/notebooks
RUN ln -s /data/ppm_rpod2 /home/user

RUN ln -s /user/cedar_scratch_andrassy /home/user/notebooks
RUN ln -s /user/cedar_scratch_andrassy /home/user
RUN ln -s /user/cedar_scratch_djstephe /home/user/notebooks
RUN ln -s /user/cedar_scratch_djstephe /home/user
RUN ln -s /user/cedar_scratch_fherwig /home/user/notebooks
RUN ln -s /user/cedar_scratch_fherwig /home/user
RUN ln -s /user/cedar_scratch_ondrea /home/user/notebooks
RUN ln -s /user/cedar_scratch_ondrea /home/user
RUN ln -s /user/cedar_scratch_pavelden /home/user/notebooks
RUN ln -s /user/cedar_scratch_pavelden /home/user

