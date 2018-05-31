FROM cyberhubs/mesahub
MAINTAINER Ed Brown

ENV DEBIAN_FRONTEND noninteractive

EXPOSE 8888

USER root

WORKDIR /root

ENV HOME /root

COPY apt_packages_dstarhub.txt /tmp/apt_packages_dstarhub.txt

RUN apt-get update && \
    apt-get -y install  curl

RUN curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash 

RUN apt-get update && \
    apt-get -y install  git-lfs


#RUN apt-get update && \
#    apt-get -y install $(cat /tmp/apt_packages_dstarhub.txt) && \
#    apt-get autoremove --yes && apt-get clean all && \
#    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN pip install --upgrade pip && \
    pip2 install --upgrade pip && \
    pip3 install --upgrade pip 

RUN echo "Installing Python packages..."

ADD ./single_pip2_requirements_dstarhub.txt /tmp/requirements/single_pip2_requirements_dstarhub.txt
ADD ./single_pip3_requirements_dstarhub.txt /tmp/requirements/single_pip3_requirements_dstarhub.txt

RUN pip2 install -r /tmp/requirements/single_pip2_requirements_dstarhub.txt
RUN pip3 install -r /tmp/requirements/single_pip3_requirements_dstarhub.txt

    
COPY jupyterhub-singleuser /usr/local/bin/jupyterhub-singleuser

COPY singleuser.sh /srv/singleuser/singleuser.sh

RUN  chmod -R 775 /home/user/notebooks


COPY .bash_aliases /home/user
COPY README_dstar.ipynb  /home/user/notebooks

RUN chown -R user:user /home/user && \
    chmod -R 775 /home/user/ && \
    chmod 775 /usr/local/bin/jupyterhub-singleuser && \
    chmod 775 /srv/singleuser/singleuser.sh

USER user
ENV HOME /home/user
WORKDIR /home/user

COPY ipython_config.py /etc/ipython/ipython_config.py

CMD ["/bin/bash", "/srv/singleuser/singleuser.sh"]

RUN ln -s /home/user/notebooks/README_dstar.ipynb /home/user
