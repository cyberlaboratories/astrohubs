# outreach-ai hub on astrohub.uvic.ca (cyberlaboratories/astrohubs, outreach-ai/hub).
# JupyterHub 6.2.0, the same version as the single-user image, pinned by digest (quay.io, 2026-11-09).
FROM quay.io/jupyterhub/jupyterhub@sha256:4700c2c99445a2e7594c57a0357c61069f11fb45cb5a23c870e903595d4d523a
RUN pip install --no-cache-dir \
      oauthenticator==17.5.0 \
      dockerspawner==14.0.0 \
      jupyterhub-idle-culler==2.0.0
