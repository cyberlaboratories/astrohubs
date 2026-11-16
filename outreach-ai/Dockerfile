# outreach-ai single-user image for astrohub (Falk, 2026-11-12/13; cyberlaboratories/astrohubs outreach-ai/image).
# The newest Jupyter stack, AI tools, and the NuGrid and MESA tooling of the old mesahub.
#   - Base: Jupyter docker-stacks scipy-notebook, built 2026-11-10 (JupyterLab 4.6.6, JupyterHub 6.2.0 single-user,
#     Python 3.13.18, Ubuntu 26.04), pinned by digest.
#   - Claude Code CLI from Anthropic's signed apt repository (stable channel). Users sign in with their own Claude
#     subscription (or an API key) by running `claude` in a terminal; the login lives in their own ~/.claude.
#   - OpenCode 2.2.5 (decision 1112-2): GitHub release binary, sha256 pinned. Updates off and OpenCode's own hosted
#     models blocked in the managed config /etc/opencode/opencode.json. Each user's ~/.config/opencode/opencode.json
#     (copied at first start, theirs to change) lists the velella spigot (default) and the Alberta (Aleph) models.
#   - Notebook Intelligence (NBI): the AI chat in JupyterLab; Claude mode runs the Claude Code CLI; set at first start to
#     the velella spigot.
#   - NuGrid: a second kernel "Python 3.11 (NuGrid, MESA)" with numpy < 2 (NuGridPy's kip_cont needs it), NuGridPy,
#     NuPyCEE (sygma, omega), JINAPyCEE, iniabu and mesa_reader; wendi-examples copied into each home at first start.
#   - MESA: the SDK's system prerequisites. MESA itself and its SDK come from the hub, mounted read-only at /opt/mesa
#     ($MESA_DIR, $MESASDK_ROOT); install_mesa.sh in the home installs other versions.
# Keys come from the hub's environment at spawn; nothing in this image is a secret.
FROM quay.io/jupyter/scipy-notebook@sha256:457ec366e08750746a80537c46111911f9658fe1c665703edb4839197dd9383c

LABEL org.opencontainers.image.title="outreach-ai" \
      org.opencontainers.image.description="Jupyter scipy-notebook 2026-11-10 + Claude Code, OpenCode, Notebook Intelligence, NuGrid and MESA tooling" \
      org.opencontainers.image.source="https://github.com/cyberlaboratories/astrohubs" \
      org.opencontainers.image.authors="arbutus-cloud-admin for Falk Herwig"

USER root

# Claude Code: Anthropic's apt repository; the build stops unless the key's fingerprint is the documented one.
# Also the MESA SDK's prerequisites (make, perl, X11 and zlib headers, a C shell, libc headers) and unzip.
ARG CLAUDE_KEY_FPR=31DDDE24DDFAB679F42D7BD2BAA929FF1A7ECACE
RUN set -eux; \
    apt-get update; \
    apt-get install -y --no-install-recommends curl gnupg ca-certificates ripgrep \
        make perl libx11-dev zlib1g-dev tcsh libc6-dev unzip; \
    install -d -m 0755 /etc/apt/keyrings; \
    curl -fsSL https://downloads.claude.ai/keys/claude-code.asc -o /etc/apt/keyrings/claude-code.asc; \
    export GNUPGHOME="$(mktemp -d)"; \
    gpg --show-keys --with-colons /etc/apt/keyrings/claude-code.asc | grep -q "^fpr:::::::::${CLAUDE_KEY_FPR}:"; \
    rm -rf "$GNUPGHOME"; unset GNUPGHOME; \
    echo "deb [signed-by=/etc/apt/keyrings/claude-code.asc] https://downloads.claude.ai/claude-code/apt/stable stable main" \
      > /etc/apt/sources.list.d/claude-code.list; \
    apt-get update; \
    apt-get install -y --no-install-recommends claude-code; \
    apt-get clean; rm -rf /var/lib/apt/lists/*; \
    claude --version

# The version in the image is the version users get; a new image brings updates.
ENV DISABLE_AUTOUPDATER=1

# OpenCode: one static binary from the GitHub release; the build stops unless its sha256 is the one pinned here
# (taken 2026-11-12 from the release's SHA256SUMS).
ARG OPENCODE_VERSION=2.2.5
ARG OPENCODE_SHA256=09f44302c10fd036912f59c81fd1aef5df45dd06f649197f9dfc9bbe2a79614e
RUN set -eux; \
    cd /tmp; \
    curl -fsSL -o opencode.tar.gz "https://github.com/sst/opencode/releases/download/v${OPENCODE_VERSION}/opencode-linux-x64.tar.gz"; \
    echo "${OPENCODE_SHA256}  opencode.tar.gz" | sha256sum -c -; \
    tar -xzf opencode.tar.gz opencode; \
    install -m 0755 opencode /usr/local/bin/opencode; \
    rm -f /tmp/opencode /tmp/opencode.tar.gz; \
    h="$(mktemp -d)"; HOME="$h" XDG_CONFIG_HOME="$h/.config" XDG_DATA_HOME="$h/.local/share" XDG_STATE_HOME="$h/.local/state" \
      XDG_CACHE_HOME="$h/.cache" opencode --version; rm -rf "$h"   # run as root: keep its files out of /home/jovyan
COPY opencode-system.json /etc/opencode/opencode.json

# NuGrid code, pinned (GitHub heads on 2026-11-13), as the old mesahub had it under /opt; wendi-examples as the skeleton
# copied into each home.
ARG NUPYCEE_SHA=1ef398beb9a26c1b53459449794ba5bfff65a3c4
ARG JINAPYCEE_SHA=5d8bf81f44ac2b46e1b56eb808fb07d52c435d45
ARG WENDI_SHA=6c703c4e32f97ffa3b5da76ed3f0e74c59abb6b9
RUN set -eux; \
    git clone -q https://github.com/NuGrid/NuPyCEE.git /opt/NuPyCEE; git -C /opt/NuPyCEE checkout -q "${NUPYCEE_SHA}"; \
    git clone -q https://github.com/becot85/JINAPyCEE.git /opt/JINAPyCEE; git -C /opt/JINAPyCEE checkout -q "${JINAPYCEE_SHA}"; \
    install -d -m 0755 /opt/outreach-ai/skel; \
    git clone -q https://github.com/NuGrid/wendi-examples.git /opt/outreach-ai/skel/wendi-examples; \
    git -C /opt/outreach-ai/skel/wendi-examples checkout -q -B master "${WENDI_SHA}"; \
    chmod -R a+rX /opt/NuPyCEE /opt/JINAPyCEE /opt/outreach-ai

# The first-start script (copies into the home only what is absent), the OpenCode preset, the shell's MESA setup.
COPY skel/ /opt/outreach-ai/skel/
RUN set -eux; \
    install -m 0755 /opt/outreach-ai/skel/10-outreach-ai-home.py /usr/local/bin/before-notebook.d/10-outreach-ai-home.py; \
    chmod 0755 /opt/outreach-ai/skel/install_mesa.sh; \
    install -m 0755 /opt/outreach-ai/skel/ai-setup /usr/local/bin/ai-setup; \
    install -m 0644 /opt/outreach-ai/skel/20-outreach-ai-keys.sh /usr/local/bin/before-notebook.d/20-outreach-ai-keys.sh; \
    install -d -m 0755 /etc/jupyter/jupyter_server_config.d; \
    install -m 0644 /opt/outreach-ai/skel/resource-usage.json /etc/jupyter/jupyter_server_config.d/outreach-ai-resource-usage.json; \
    cat /opt/outreach-ai/skel/bash.bashrc.outreach-ai >> /etc/bash.bashrc; \
    install -d -m 0755 /opt/outreach-ai/skel/dot; \
    for f in .bashrc .profile .bash_logout; do [ -f "/home/${NB_USER}/$f" ] && cp "/home/${NB_USER}/$f" /opt/outreach-ai/skel/dot/; done; \
    chmod -R a+rX /opt/outreach-ai; \
    "${CONDA_DIR}/bin/python" -c "import json, os; p = '${CONDA_DIR}/share/jupyter/lab/settings/overrides.json'; d = json.load(open(p)) if os.path.exists(p) else {}; d.setdefault('@jupyterlab/docmanager-extension:plugin', {}).setdefault('defaultViewers', {})['markdown'] = 'Markdown Preview'; d.setdefault('@jupyterlab/apputils-extension:notification', {}).update({'fetchNews': 'false', 'checkForUpdates': False}); os.makedirs(os.path.dirname(p), exist_ok=True); json.dump(d, open(p, 'w'), indent=2); print('lab overrides:', json.dumps(d))"; \
    fix-permissions "${CONDA_DIR}/share/jupyter/lab/settings"   # .md opens as preview; no "Jupyter news" prompt and no update check (Falk); root wrote it, so hand it back

USER ${NB_UID}

# Node.js for MCP servers started with npx (NBI); Notebook Intelligence pinned; mesa_reader in the main Python too.
ARG NBI_VERSION=6.3.0
RUN set -eux; \
    mamba install -y -c conda-forge nodejs; \
    mamba clean -afy; \
    mamba install -y -c conda-forge gh; \
    mamba clean -afy; \
    pip install --no-cache-dir "notebook-intelligence==${NBI_VERSION}" "mesa-reader==0.4.0" "iniabu==1.2.1" "ads==0.13.0" \
        "astropy==8.0.2" "astroquery==0.4.11" \
        "jupyter-resource-usage==1.3.0" \
        "jupyterlab-lsp==5.3.0" "python-lsp-server==1.15.0" "pyflakes==4.1.0"; \
    conda config --system --set changeps1 False; \
    jupyter labextension disable --level=sys_prefix @jupyterlab/extensionmanager-extension; \
    fix-permissions "${CONDA_DIR}"; fix-permissions "/home/${NB_USER}"

# The NuGrid kernel: Python 3.11 with numpy < 2. NuGridPy at the GitHub head (0.8.0); iniabu 1.1.2, the last release
# without numpy >= 2 (1.2.x requires it); NuPyCEE and JINAPyCEE on the kernel's PYTHONPATH from /opt, as in the old
# mesahub (they carry their data tables in the checkout).
ARG NUGRIDPY_SHA=7202f6474119393de82b5eb0670bb08cc321e4de
RUN set -eux; \
    mamba create -y -n nugrid -c conda-forge python=3.11 "numpy<2" scipy matplotlib h5py pandas seaborn ipykernel \
        future xlrd ipywidgets; \
    "${CONDA_DIR}/envs/nugrid/bin/pip" install --no-cache-dir \
        "git+https://github.com/NuGrid/NuGridPy.git@${NUGRIDPY_SHA}" "iniabu==1.1.2" "mesa-reader==0.4.0" "numpy<2"; \
    "${CONDA_DIR}/envs/nugrid/bin/python" -c "import numpy, sys; v = numpy.__version__; print('nugrid numpy', v); sys.exit(int(v.split('.')[0]) >= 2)"; \
    "${CONDA_DIR}/envs/nugrid/bin/python" -m ipykernel install --prefix "${CONDA_DIR}" --name nugrid \
        --display-name "Python 3.11 (NuGrid, MESA)"; \
    "${CONDA_DIR}/envs/nugrid/bin/python" -c "import json, sys; p = sys.argv[1]; k = json.load(open(p)); k['env'] = {'PYTHONPATH': '/opt/NuPyCEE:/opt/JINAPyCEE:/opt'}; json.dump(k, open(p, 'w'), indent=1)" \
        "${CONDA_DIR}/share/jupyter/kernels/nugrid/kernel.json"; \
    mamba clean -afy; \
    fix-permissions "${CONDA_DIR}"; fix-permissions "/home/${NB_USER}"

# The hub runs the user as "user" with home /home/user and mounts <name>_config INSIDE that home. Without a /home/user in
# the image, Docker would create it for the mount, owned by root and empty, and the user could not write their home.
# So the image carries a ready /home/user (a copy of jovyan's), owned by the user; docker-stacks then just points the
# renamed user at it.
USER root
RUN cp -a "/home/${NB_USER}" /home/user && chown -R "${NB_UID}:${NB_GID}" /home/user
USER ${NB_UID}

WORKDIR /home/${NB_USER}

# Notebook Intelligence 6.3.0: the start-folder picker becomes a text box (Enter in JupyterLab's FileDialog closed it mid-rename: "File browser is disposed"; card 1113-1). Fails the build if NBI changes.
USER root
COPY patches/nbi-folder-picker.py /opt/outreach-ai/patches/nbi-folder-picker.py
RUN python3 /opt/outreach-ai/patches/nbi-folder-picker.py
# Notebook Intelligence 6.3.0: the Claude Code launcher tile lists earlier sessions even when the chat panel's Claude mode is off (card 1113-6). Fails the build if NBI changes.
COPY patches/nbi-claude-sessions.py /opt/outreach-ai/patches/nbi-claude-sessions.py
RUN python3 /opt/outreach-ai/patches/nbi-claude-sessions.py
USER ${NB_UID}
WORKDIR /home/${NB_USER}
