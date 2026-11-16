# Outreach-AI: NuGrid, MESA and AI assistants on Astrohub

Outreach-AI is a JupyterHub on [Astrohub](https://astrohub.uvic.ca), run by the [Computational Stellar Astrophysics](https://www.ppmstar.org) group in the [Department of Physics and Astronomy](https://www.uvic.ca/science/physics) at the [University of Victoria](https://www.uvic.ca). It has the current Jupyter stack (JupyterLab 4.6, Python 3.13), three AI assistants, the NuGrid Python tools and the MESA stellar evolution code. We are trying it out with a few students.

### Index

* [Access and file storage](#storage)
* [AI assistants](#ai)
* [NuGrid](#nugrid)
* [MESA](#mesa)
* [Who runs this hub](#who)

<a id='storage'></a>
### Access and file storage

* You log in with GitHub. Access is by invitation.
* As on the other Astrohub hubs, the data in your home directory lives inside your container. It is lost when you or the administrator stop your server (via `File` → `Hub Control Panel`), but not when you just log out or close the browser. A server idle for one week is stopped.
* `yourgithubname_config` in your home directory is your private persistent folder, the same one you see on the Outreach hub. Keep keys and config files there.
* `yourgithubname_config/outreach-ai-ini` keeps your logins and settings for this hub: the Claude login, ssh keys, git settings, OpenCode settings and sessions, the chat settings in JupyterLab, and the keys you store with `ai-setup`. Every server start links them into your new home directory, so you enter them only once. A file `bashrc` that you put there is read by every new terminal.
* Keep your work in `/user/scratch14_outreach` (see below) or in a git repository you push often.
* The file browser opens in your home directory, `/home/user`; your unix name is `user`, as on the other Astrohub hubs. Its top level is `/`, so that `/data` and `/user` can be reached; your home directory also has the links `data` and `user`.

Directory | Content
----------|--------
`/data` | Read-only data repositories: `ASDR` (Astrophysics Simulation Data Repository, with the NuGrid stellar evolution and yield data sets, e.g. Ritter et al. 2018, 2018MNRAS.480..538R) and `nugrid_data` (the NuGrid data repository, also on Globus as `astrohub#nugrid`)
`/user` | The shared space of the Outreach hub, read-only here (pre-installed MESA versions and more)
`/user/scratch14_outreach` | Writable shared scratch space, the same as on the Outreach hub
`/opt/mesa` | MESA and the MESA SDK, read-only (see below)

* Everyone on this hub and on the Outreach hub works in `/user/scratch14_outreach` as the same unix user. This is on purpose, so that people can collaborate without barriers. Create a directory with your name and keep to it, do not change other people's files, and keep a copy elsewhere of anything you cannot lose.

* The machine has 4 CPUs and 29 GB of memory, shared by everyone on Astrohub. Your server has no fixed limit, so please stop large runs you no longer need.

<a id='ai'></a>
### AI assistants

#### Checking your setup

In a terminal, `ai-setup` checks your account and says what, if anything, to do: the Velella and Aleph models, OpenCode and the chat in JupyterLab, your Claude login, git and ssh, your own keys, and MESA. `ai-setup fix` walks you through what is missing: your git name and email, an ssh key, and the Claude login. What you set up is kept in `yourgithubname_config/outreach-ai-ini`, so it is there at every later login.

#### Your own keys

`ai-setup key add` stores a key of your own: a GitHub token (`GH_TOKEN`, used by the `gh` command), a NASA ADS token (`ADS_API_TOKEN`, used by the `ads` Python package), an OpenAI, Anthropic, Gemini or OpenRouter key (OpenCode uses these by itself), or any other by name. It asks for the value without showing it and keeps it in `~/.config/outreach-ai/keys.env`, readable only by you. New terminals have your keys at once, notebooks after your next server start. `ai-setup key list` shows the names, `ai-setup key remove NAME` deletes one. `ai-setup model add` adds your own OpenAI-compatible model provider to OpenCode. An Anthropic API key makes Claude Code bill that key instead of your subscription.

#### Installing software

Anything you install (`pip install --user`, a JupyterLab extension) is lost when your server stops, which is why the JupyterLab extension manager is switched off here. For a Python environment of your own that lasts, create it in `/user/scratch14_outreach`; ask Falk if a package or extension should be in the image for everyone.

#### Claude Code

In a terminal, run `claude` and log in with your own Claude account. A Pro or Max subscription, or an Anthropic API key, is needed; free accounts cannot use Claude Code. The login is kept in `yourgithubname_config/outreach-ai-ini/claude`, linked as `~/.claude`. The chat in JupyterLab can use the same login: open its settings (the gear in the chat panel) and turn on Claude mode.

#### OpenCode

In a terminal, run `opencode`. It is set up with these open-weight models; switch with `/models`:

Model | Provider | Context
------|----------|--------
`velella/qwen3.6-35b-a3b` (default) | Velella spigot | 64k tokens
`aleph/qwen38-27b` | Aleph, University of Alberta | 512k tokens
`aleph/qwen35-122b` | Aleph | 128k tokens
`aleph/gemma-4-26b-a4b` | Aleph | 128k tokens
`aleph/gpt-oss-120b` | Aleph | 128k tokens
`aleph/gpt-oss-20b` | Aleph | 128k tokens

Aleph serves more models than these, for example `astrosage`. To add one, or your own provider, edit `~/.config/opencode/opencode.json`; it is yours to change. OpenCode's own hosted models are turned off on this hub.

#### The chat in JupyterLab

The chat (Notebook Intelligence) is set to the Velella model. To use an Aleph model instead, open its settings and change only the model id, for example to `qwen38-27b` or `gpt-oss-120b`; keep the base URL.

#### How the models are reached

OpenCode and the chat reach the Velella and Aleph models through the hub's model gateway, `http://outreach-ai-gateway:8080/v1`, which holds the hub's keys. Your tools sign in to it with your own hub login (the `JUPYTERHUB_API_TOKEN` your server has), so there are no shared keys in your account, and the models work only from this hub. The Velella spigot serves course students first, so it can be slower when a class is working.

#### Where the assistants work

The assistants read and change files. Start them in a folder of your own, for example `/user/scratch14_outreach/yourgithubname`, and take care: others share `/user/scratch14_outreach`.

<a id='nugrid'></a>
### NuGrid

* The kernel "Python 3.11 (NuGrid, MESA)" has [NuGridPy](https://github.com/NuGrid/NuGridPy), [NuPyCEE](https://github.com/NuGrid/NuPyCEE) (`sygma`, `omega`), [JINAPyCEE](https://github.com/becot85/JINAPyCEE), `iniabu` and `mesa_reader`, with NumPy below 2, as NuGridPy needs. Use it for the NuGrid notebooks. The default Python 3 kernel has the scientific Python stack and `mesa_reader`.
* `wendi-examples` in your home directory is a copy of [NuGrid/wendi-examples](https://github.com/NuGrid/wendi-examples); run `git pull` in it for the latest version. `Star_explore.ipynb` is a good place to start.
* Some examples refer to data paths that exist only on other hubs.

<a id='mesa'></a>
### MESA

[MESA](https://docs.mesastar.org) r25.12.1 is installed at `$MESA_DIR`, with its SDK at `$MESASDK_ROOT`, and every terminal loads the SDK. `$MESA_DIR` is read-only, and your home directory is lost when the server stops, so run MESA in a folder of your own in `/user/scratch14_outreach`.

Your first run, in a terminal (replace `yourgithubname`):

```
mkdir -p /user/scratch14_outreach/yourgithubname && cd /user/scratch14_outreach/yourgithubname
cp -r $MESA_DIR/star/work mesa-test && cd mesa-test
export OMP_NUM_THREADS=4
./mk && ./rn
```

`./mk` compiles your run in about a minute, and `./rn` evolves the default model. The output goes to `LOGS` and `photos` in that folder and is kept when your server stops. To plot it, use the kernel "Python 3.11 (NuGrid, MESA)", for example with `mesa_reader` or NuGridPy's `mesa` module.

* `OMP_NUM_THREADS` is 2 unless you set it; 4 uses the whole machine, so lower it again when others are working.
* To run a test-suite case, copy it out the same way (for example `cp -r $MESA_DIR/star/test_suite/1M_pre_ms_to_wd .`, a Sun-like star to a white dwarf) and follow the MESA documentation on running a test case outside the test suite.
* For another MESA version, `~/install_mesa.sh` downloads a version and its matching SDK into `/user/scratch14_outreach/yourgithubname/mesa` and compiles it, for example `~/install_mesa.sh 24.08.1`. Each version takes several GB and up to an hour. `~/install_mesa.sh --list` shows the versions it knows. Only 25.12.1 has been built on this hub so far. Recent versions should work; older ones need older SDKs, which may not run on this hub's newer Linux, and the oldest also need fixes to the MESA source. The versions in `/user/mesa` were built for the other hubs.

<a id='who'></a>
### Who runs this hub

Falk Herwig and the CSA group at UVic. The code and this README are in [cyberlaboratories/astrohubs](https://github.com/cyberlaboratories/astrohubs), folder `outreach-ai`. Questions and problems go to Falk.
