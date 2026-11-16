# outreach-ai: the user's own keys (stored with `ai-setup key add`) into the Jupyter server's environment, so notebooks
# and their kernels have them. docker-stacks SOURCES this from before-notebook.d, as root when the hub starts the
# container as root. So the file is READ line by line, never run: only lines NAME=value with an upper-case NAME are taken,
# and names the hub sets are skipped.
# The persistent place for logins and settings, made by the first-start script in <name>_config (the home itself is
# temporary, as on the other Astrohub hubs). Claude Code keeps its whole config in ~/.claude (a link into it), git its
# settings in a file there.
_oa_home="/home/${NB_USER:-jovyan}"
if [ -n "${JUPYTERHUB_USER:-}" ] && [ -d "$_oa_home/${JUPYTERHUB_USER}_config/outreach-ai-ini" ]; then
  export OUTREACH_AI_INI="$_oa_home/${JUPYTERHUB_USER}_config/outreach-ai-ini"
  export CLAUDE_CONFIG_DIR="$_oa_home/.claude"
  export GIT_CONFIG_GLOBAL="$OUTREACH_AI_INI/gitconfig"
fi
unset _oa_home
_oa_keys="/home/${NB_USER:-jovyan}/.config/outreach-ai/keys.env"
if [ -r "$_oa_keys" ]; then
  while IFS= read -r _oa_line || [ -n "$_oa_line" ]; do
    _oa_name=${_oa_line%%=*}
    case "$_oa_name" in
      JUPYTERHUB_API_TOKEN|OUTREACH_AI_GATEWAY|OUTREACH_AI_INI|CLAUDE_CONFIG_DIR|GIT_CONFIG_GLOBAL|VELA_SPIGOT_API_KEY|ALEPH_API_KEY|PATH|HOME|USER|SHELL|MESA_DIR|MESASDK_ROOT|OMP_NUM_THREADS|NB_USER|NB_UID|NB_GID) continue ;;
    esac
    [[ "$_oa_name" =~ ^[A-Z_][A-Z0-9_]*$ ]] && [ "$_oa_name" != "$_oa_line" ] && export "$_oa_name=${_oa_line#*=}"
  done < "$_oa_keys"
fi
unset _oa_keys _oa_line _oa_name
