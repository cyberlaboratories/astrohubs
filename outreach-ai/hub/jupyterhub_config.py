"""outreach-ai hub on astrohub.uvic.ca (cyberlaboratories/astrohubs, outreach-ai/hub).

JupyterHub 6.2 behind the Astrohub nginx at /outreach-ai/. GitHub login, admitted by an allow list that is re-read at
every login (so adding a student needs no restart). DockerSpawner, one container per user, no memory or CPU limit (each
user may use the whole machine while on), a 10 GB cap on the container's writable layer only, idle servers stopped after
one week. Secrets come from the environment (the host's .env, mode 0600); nothing secret is in this file.
"""
import os
import sys

from oauthenticator.github import GitHubOAuthenticator

c = get_config()  # noqa: F821  (provided by JupyterHub)

WLIST = "/srv/jupyterhub/access/wlist"    # allow list: GitHub user names, one per line; '#' starts a comment
BLIST = "/srv/jupyterhub/access/blist"    # block list, same format; wins over the allow list (as on the other hubs)
USERS = "/jupyterhub_users"               # /mnt/jupyterhub_users on the host (CephFS)


def names(path):
    try:
        with open(path) as fh:
            return {ln.split("#")[0].strip().lower() for ln in fh if ln.split("#")[0].strip()}
    except FileNotFoundError:
        return set()


class OutreachAIAuthenticator(GitHubOAuthenticator):
    """GitHub login, admitted by the access files re-read at the moment of login, as on the other Astrohub hubs:
    a name in wlist and not in blist, or an admin."""

    async def check_allowed(self, username, auth_model):
        if username in self.admin_users:
            return True
        name = username.lower()
        return name in names(WLIST) and name not in names(BLIST)


# --- login --------------------------------------------------------------------------------------------------------
c.JupyterHub.authenticator_class = OutreachAIAuthenticator
c.GitHubOAuthenticator.client_id = os.environ["GITHUB_CLIENT_ID"]
c.GitHubOAuthenticator.client_secret = os.environ["GITHUB_CLIENT_SECRET"]
c.GitHubOAuthenticator.oauth_callback_url = "https://astrohub.uvic.ca/outreach-ai/hub/oauth_callback"
c.Authenticator.admin_users = {n.strip().lower() for n in os.environ.get("OUTREACH_AI_ADMINS", "").split(",") if n.strip()}
c.Authenticator.allow_existing_users = False      # the allow list decides, also for users seen before

# --- the hub itself -----------------------------------------------------------------------------------------------
c.JupyterHub.base_url = "/outreach-ai/"
c.JupyterHub.bind_url = "http://0.0.0.0:8000/outreach-ai/"   # compose maps 127.0.0.1:10443 to this
c.JupyterHub.hub_ip = "0.0.0.0"                               # the hub API, for the users' containers
c.JupyterHub.hub_connect_ip = "multiuser_outreach_ai"         # the hub's name on outreach-ai-net
c.JupyterHub.cookie_secret_file = "/srv/jupyterhub/state/jupyterhub_cookie_secret"
c.JupyterHub.cleanup_servers = False   # a hub restart leaves users' servers running (default True stopped them, 11-13)
c.JupyterHub.db_url = "sqlite:////srv/jupyterhub/state/jupyterhub.sqlite"


# --- the users' containers ----------------------------------------------------------------------------------------
def make_home(spawner):
    """Create the user's persistent <name>_config folder on first login, as the outreach hub makes it (same folder)."""
    top = os.path.join(USERS, spawner.user.name)
    if not os.path.isdir(top):
        os.mkdir(top, 0o755)
        os.chown(top, 1000, 1000)


c.Spawner.pre_spawn_hook = make_home
c.JupyterHub.spawner_class = "dockerspawner.DockerSpawner"
c.DockerSpawner.image = os.environ.get("OUTREACH_AI_IMAGE", "local/outreach:ai-2026-11-13")
c.DockerSpawner.network_name = "outreach-ai-net"
c.DockerSpawner.use_internal_ip = True
c.DockerSpawner.remove = True
c.DockerSpawner.name_template = "outreach-ai-{username}"
# The file browser's top is / so that /data and /user can be reached (Jupyter does not follow links out of its top
# folder); it opens at home, and the home has data and user links as on the other hubs (Falk, 2026-11-13).
c.DockerSpawner.notebook_dir = "/"
c.Spawner.args = ["--ServerApp.preferred_dir=/home/user"]
# The unix user is "user" with home /home/user, as on the old Astrohub images, and the prompt reads user@outreach-ai
# (Falk, 2026-11-13). docker-stacks renames its default user jovyan when started as root with NB_USER, then drops to it.
c.DockerSpawner.extra_create_kwargs = {"user": "root", "hostname": "outreach-ai"}
# As on the other Astrohub hubs (Falk, 2026-11-13): the home lives in the container and is lost when the server stops;
# <name>_config is the persistent folder, the same one the outreach hub mounts. The image keeps the user's logins and
# settings in <name>_config/outreach-ai-ini and links them into each new home.
c.DockerSpawner.volumes = {
    "/mnt/jupyterhub_users/{username}": "/home/user/{username}_config",
    "/mnt/outreach/data": {"bind": "/data", "mode": "ro"},
    "/mnt/outreach/user": {"bind": "/user", "mode": "ro"},
    # the shared user scratch space, writable as on the outreach hub; the rest of /user stays read-only
    "/mnt/outreach/user/scratch14_outreach": {"bind": "/user/scratch14_outreach", "mode": "rw"},
    "/mnt/jupyter_config/outreach-ai": {"bind": "/jupyterhub_config", "mode": "ro"},
    "/mnt/outreach/mesa-ai": {"bind": "/opt/mesa", "mode": "ro"},
}
# No memory or CPU limit (Falk, 2026-11-13). The one cap: the writable layer, because Docker's disk is shared by every
# hub on astrohub; homes and MESA are on CephFS and do not count against it.
c.DockerSpawner.extra_host_config = {"storage_opt": {"size": "10G"}}
c.DockerSpawner.environment = {
    "NB_USER": "user", "NB_UID": "1000", "NB_GID": "100",
    # No model key enters a user's server (Falk, 2026-11-13: "can't we hide them?"): the gateway container holds them and
    # users' tools reach it with their own JUPYTERHUB_API_TOKEN, which the gateway checks with the hub.
    "OUTREACH_AI_GATEWAY": "http://outreach-ai-gateway:8080/v1",
    "MESA_DIR": "/opt/mesa/mesa-25.12.1",          # MESA zips from r24.08.1 on unpack without the "r"
    "MESASDK_ROOT": "/opt/mesa/sdk-25.12.1/mesasdk",   # the layout install_mesa.sh makes
    "OMP_NUM_THREADS": "2",
}
c.Spawner.default_url = "/lab"
c.Spawner.start_timeout = 300
c.Spawner.http_timeout = 180

# --- stop servers idle for one week (Falk, 2026-11-13) -------------------------------------------------------------
c.JupyterHub.services = [{
    "name": "idle-culler",
    "command": [sys.executable, "-m", "jupyterhub_idle_culler", "--timeout=604800", "--cull-every=3600"],
}]
c.JupyterHub.load_roles = [{
    "name": "idle-culler",
    "scopes": ["list:users", "read:users:activity", "read:servers", "delete:servers"],
    "services": ["idle-culler"],
}]
