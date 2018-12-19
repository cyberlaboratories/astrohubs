FROM local/basichub
MAINTAINER Falk Herwig

ENV DEBIAN_FRONTEND noninteractive

EXPOSE 8888

USER root

WORKDIR /root

ENV HOME /root

RUN pip install --upgrade pip && \
    pip2 install --upgrade pip && \
    pip3 install --upgrade pip 

RUN echo "Installing Python packages ..."

ADD ./single_pip2_requirements_ppmstar.txt /tmp/requirements/single_pip2.txt
ADD ./single_pip3_requirements_ppmstar.txt /tmp/requirements/single_pip3.txt

RUN pip2 install -r /tmp/requirements/single_pip2.txt  
RUN pip3 install -r /tmp/requirements/single_pip3.txt 

COPY jupyterhub-singleuser /usr/local/bin/jupyterhub-singleuser
COPY singleuser.sh /srv/singleuser/singleuser.sh
#RUN  chmod -R 775 /home/user/notebooks

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


COPY .bash_aliases /home/user
RUN echo restart here
COPY README_ppmstar.ipynb  /home/user/

RUN echo "Loading PPMnotebooks git repo ... "
RUN git clone https://github.com/PPMstar/PPMnotebooks.git /home/user/PPMnotebooks

RUN chown -R user:user /home/user && \
    chmod -R 775 /home/user/ && \
    chmod 775 /usr/local/bin/jupyterhub-singleuser && \
    chmod 775 /srv/singleuser/singleuser.sh

USER user
ENV HOME /home/user
WORKDIR /home/user

CMD ["/bin/bash", "/srv/singleuser/singleuser.sh"]


