#!/opt/conda/bin/python
"""outreach-ai: set up the user's home at server start (docker-stacks runs this from /usr/local/bin/before-notebook.d).

Copies in only what is ABSENT, so a user's own changes are never overwritten:
  - links from the (temporary) home into <name>_config/outreach-ai-ini, the persistent place for logins and settings,
  - shell dot files, wendi-examples, links to the hub README, Star_explore.ipynb, /data and /user, install_mesa.sh,
  - a MESA starter work folder (when the hub mounts MESA at $MESA_DIR),
  - OpenCode settings (velella by default, the Alberta models listed), through the hub's model gateway,
  - the Notebook Intelligence chat settings, through the gateway with the user's hub token (refreshed at each start).
  The shared model keys live only in the gateway, never in a user's server.
Never fails the server start: every step is wrapped, and problems are printed to the container log.
"""
import json
import os
import shutil
import sys

if os.geteuid() == 0:
    # The hub starts the container as root with NB_USER=user (the name the old Astrohub images used); docker-stacks
    # renames jovyan and then runs these hooks still as root. Set the home up as that user, so nothing in it is root's.
    os.execvp("sudo", ["sudo", "-E", "-H", "-u", os.environ.get("NB_USER", "jovyan"),
                       sys.executable, os.path.abspath(__file__)])

HOME = os.path.expanduser("~")
SKEL = "/opt/outreach-ai/skel"
NAME = os.environ.get("JUPYTERHUB_USER", "")
CFG = os.path.join(HOME, f"{NAME}_config") if NAME else ""      # persistent, the same folder the outreach hub mounts
INI = os.path.join(CFG, "outreach-ai-ini") if CFG else ""        # logins and settings kept across servers
# home path -> folder inside outreach-ai-ini. The home itself is temporary (as on the other Astrohub hubs); these are
# linked into each new home so a login or a key entered once stays.
KEPT = ((".claude", "claude"), (".ssh", "ssh"), (".config/opencode", "opencode"),
        (".local/share/opencode", "opencode-data"), (".jupyter/nbi", "nbi"), (".config/outreach-ai", "outreach-ai"))
GATEWAY = os.environ.get("OUTREACH_AI_GATEWAY", "http://outreach-ai-gateway:8080/v1")   # holds the shared keys
OLD_URLS = {"https://tokens.velella.ca/v1", "https://inference.vulcan.alliancecan.ca/v1"}  # direct, before the gateway
DEFAULT_MODEL = "qwen3.6-35b-a3b"


def step(name, fn):
    try:
        fn()
    except Exception as e:  # never block the server
        print(f"outreach-ai home: {name}: {type(e).__name__}: {e}", file=sys.stderr)


def copy_if_absent(src, dst):
    if os.path.lexists(dst) or not os.path.exists(src):
        return
    if os.path.isdir(src):
        shutil.copytree(src, dst, symlinks=True)
    else:
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copy2(src, dst)


def link_if_absent(target, dst):
    if not os.path.lexists(dst) and os.path.exists(target):
        os.symlink(target, dst)


def persist():
    """Link the kept folders from <name>_config/outreach-ai-ini into the home (made on first use)."""
    if not CFG or not os.path.isdir(CFG):
        return                                   # no persistent folder mounted (a test container): nothing to keep
    os.makedirs(INI, exist_ok=True)
    for home_rel, ini_rel in KEPT:
        target = os.path.join(INI, ini_rel)
        os.makedirs(target, exist_ok=True)
        link = os.path.join(HOME, home_rel)
        if os.path.islink(link):
            continue
        if os.path.isdir(link):                  # a real folder from the image: keep anything new, then link
            for n in os.listdir(link):
                if not os.path.exists(os.path.join(target, n)):
                    shutil.move(os.path.join(link, n), os.path.join(target, n))
            shutil.rmtree(link)
        os.makedirs(os.path.dirname(link), exist_ok=True)
        os.symlink(target, link)
    os.chmod(os.path.join(INI, "ssh"), 0o700)
    os.chmod(os.path.join(INI, "outreach-ai"), 0o700)
    readme = os.path.join(INI, "README.md")
    if not os.path.exists(readme):
        with open(readme, "w") as fh:
            fh.write("# outreach-ai-ini\n\nYour outreach-ai logins and settings, kept here because the home directory is "
                     "temporary. Each server start links these folders into the new home: claude (Claude Code login), ssh, "
                     "opencode and opencode-data, nbi (the JupyterLab chat), outreach-ai (keys stored with ai-setup). "
                     "gitconfig is your git settings; a file named bashrc here is read by every new terminal.\n")


def dotfiles():
    for f in (".bashrc", ".profile", ".bash_logout"):
        copy_if_absent(os.path.join(SKEL, "dot", f), os.path.join(HOME, f))


def examples():
    copy_if_absent(os.path.join(SKEL, "wendi-examples"), os.path.join(HOME, "wendi-examples"))
    link_if_absent("/jupyterhub_config/README.md", os.path.join(HOME, "README.md"))
    link_if_absent("/data", os.path.join(HOME, "data"))     # as on the other Astrohub hubs
    link_if_absent("/user", os.path.join(HOME, "user"))
    link_if_absent("wendi-examples/Stellar evolution and nucleosynthesis data/Star_explore.ipynb",
                   os.path.join(HOME, "Star_explore.ipynb"))
    copy_if_absent(os.path.join(SKEL, "install_mesa.sh"), os.path.join(HOME, "install_mesa.sh"))


def mesa_work():
    mesa = os.environ.get("MESA_DIR", "")
    if mesa and os.path.isdir(os.path.join(mesa, "star", "work")):
        copy_if_absent(os.path.join(mesa, "star", "work"), os.path.join(HOME, "mesa-work"))


def opencode():
    d = os.path.join(HOME, ".config", "opencode")
    if not any(os.path.exists(os.path.join(d, n)) for n in ("config.json", "opencode.json", "opencode.jsonc")):
        copy_if_absent(os.path.join(SKEL, "opencode.json"), os.path.join(d, "opencode.json"))
        return
    # Settings made before the gateway (2026-11-13) still point straight at Velella or Aleph with a shared key that the
    # server no longer has: move just those two providers to the gateway and the user's hub token. Nothing else changes.
    p = os.path.join(d, "opencode.json")
    if not os.path.exists(p):
        return
    with open(p) as fh:
        cfg = json.load(fh)
    changed = False
    for prov in (cfg.get("providers") or {}).values():
        s = prov.get("settings") or {}
        if s.get("baseURL") in OLD_URLS:
            s["baseURL"] = GATEWAY
            prov["env"] = ["JUPYTERHUB_API_TOKEN"]
            changed = True
    if changed:
        with open(p + ".tmp", "w") as fh:
            json.dump(cfg, fh, indent=2)
        os.replace(p + ".tmp", p)


def nbi():
    # The chat reaches the models through the hub's gateway with the user's hub token. That token changes at every
    # server start and NBI keeps its key in its own config (no env fallback for this provider), so the key is refreshed
    # here at each start whenever the config points at the gateway; configs from before the gateway are moved to it.
    key = os.environ.get("JUPYTERHUB_API_TOKEN", "")
    path = os.path.join(HOME, ".jupyter", "nbi", "config.json")

    def model(kind):
        return {"provider": "openai-compatible", "model": f"openai-compatible-{kind}-model",
                "properties": [{"id": "base_url", "value": GATEWAY}, {"id": "api_key", "value": key},
                               {"id": "model_id", "value": DEFAULT_MODEL}]}

    if not os.path.exists(path):
        if not key:
            return
        os.makedirs(os.path.dirname(path), exist_ok=True)
        cfg = {"chat_model": model("chat"), "inline_completion_model": model("inline-completion")}
    else:
        with open(path) as fh:
            cfg = json.load(fh)
        changed = False
        for k in ("chat_model", "inline_completion_model"):
            props = (cfg.get(k) or {}).get("properties") or []
            vals = {p.get("id"): p for p in props if isinstance(p, dict)}
            if vals.get("base_url", {}).get("value") in OLD_URLS:
                vals["base_url"]["value"] = GATEWAY
                changed = True
            if key and vals.get("base_url", {}).get("value") == GATEWAY and "api_key" in vals \
                    and vals["api_key"].get("value") != key:
                vals["api_key"]["value"] = key
                changed = True
        if not changed:
            return
    tmp = path + ".tmp"
    with open(os.open(tmp, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600), "w") as fh:
        json.dump(cfg, fh, indent=2)
    os.replace(tmp, path)


for name, fn in (("dot files", dotfiles), ("kept settings", persist), ("examples", examples),
                 ("OpenCode", opencode), ("NBI", nbi)):   # MESA runs go in /user/scratch14_outreach (README), not the home
    step(name, fn)
