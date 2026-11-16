#!/bin/bash
# install_mesa.sh: download, check and compile a MESA version together with its matching MESA SDK.
# For the outreach-ai hub on Astrohub (cyberlaboratories/astrohubs, outreach-ai); successor of mesahub/install_mesa.sh.
#
#   ~/install_mesa.sh --list                 the versions this script knows
#   ~/install_mesa.sh 24.08.1                install into /user/scratch14_outreach/<you>/mesa (kept when the server stops)
#   ~/install_mesa.sh 24.08.1 /some/dir      install somewhere else
#
# Layout: <dir>/mesa-<version>/ for MESA and <dir>/sdk-<sdk version>/ for its SDK. Downloads are checked against Zenodo's
# md5 before anything is unpacked. Each version needs a few GB and up to an hour to compile; the log goes to
# <MESA_DIR>/install.log. At the end it prints the three lines that select this version in your shell.
# Tested end to end: 25.12.1. The others follow the same steps but have not been compiled on this hub yet.
set -euo pipefail

# version | zenodo record | zip file | md5 | folder in the zip | SDK version | SDK record | SDK file | SDK md5 | SDK folder
TABLE="
26.04.1|19722306|mesa-26.04.1.zip|994fd28f38da92efc80e38ca11d8317b|mesa-26.04.1|26.6.1|20598423|mesasdk-x86_64-linux-26.6.1.tar.gz|c5ef95aafe07848988959fd45cfe8baf|mesasdk-26.6.1
25.12.1|18023257|mesa-25.12.1.zip|f5c53511f8ac0114bb0025558fde8e57|mesa-25.12.1|25.12.1|18163992|mesasdk-x86_64-linux-25.12.1.tar.gz|fd6de5df46c47b63fb1e2b18ed3dd884|mesasdk
24.08.1|13353788|mesa-24.08.1.zip|75418c76ad48783cdb8a2cfdc8be020f|mesa-24.08.1|24.7.1|13768913|mesasdk-x86_64-linux-24.7.1.tar.gz|77939513276648cda868e9a0c9cdd025|mesasdk
23.05.1|7983526|mesa-r23.05.1.zip|77d598cc6db7e1714c0c2ecab61973ad|mesa-r23.05.1|23.7.3|10624843|mesasdk-x86_64-linux-23.7.3.tar.gz|82ba6f4da21c8afe8e50bb0bce1b9435|mesasdk
22.05.1|6547951|mesa-r22.05.1.zip|628be12fa7f4ab31de43646c5de591c9|mesa-r22.05.1|22.6.1|7457681|mesasdk-x86_64-linux-22.6.1.tar.gz|87d4fc91ce9eea631739588206ec4fff|mesasdk
21.12.1|5798242|mesa-r21.12.1.zip|fc2b3692b98391be83ab8251ee354d42|mesa-r21.12.1|21.4.1|5802444|mesasdk-x86_64-linux-21.4.1.tar.gz|e548fa0e107f08bd7ce5ce298ec96d61|mesasdk
15140|4311514|mesa-r15140.zip|fcd8ef3dd72ff25890ef2eb8201e9f8f|mesa-r15140|20.12.1|4587206|mesasdk-x86_64-linux-20.12.1.tar.gz|e1d65193f2c78557f50c02289ea08ada|mesasdk
12778|3698354|mesa-r12778.zip|314c0ac396a9540d88bfa9b3e63574d1|mesa-r12778|20.3.2|4587200|mesasdk-x86_64-linux-20.3.2.tar.gz|e5027b3fcda99d7800af1f7c1bcc80da|mesasdk
12115|3473377|mesa-r12115.zip|e8baf602d9b9cc1c301a404d4acc2215|mesa-r12115|20190830|3560834|mesasdk-x86_64-linux-20190830.tar.gz|8d618120df849d6a04c479544965c137|mesasdk
10398|2630891|mesa-r10398.zip|6441be1849e8dfefe6d0aa054b9385b0|mesa-r10398|20180822|2603170|mesasdk-x86_64-linux-20180822.tar.gz|d4312b3eed6f66efc5d5f1f64789f6a7|mesasdk
"

usage() { sed -n '2,13p' "$0" | sed 's/^# \{0,1\}//'; }
[ $# -ge 1 ] || { usage; exit 1; }
if [ "$1" = --list ] || [ "$1" = -l ]; then
  echo "MESA version  (SDK)"; echo "$TABLE" | awk -F'|' 'NF {printf "  %-10s  (%s)\n", $1, $6}'; exit 0
fi
want=${1#r}                       # accept r15140 and 15140, r24.08.1 and 24.08.1
# Default: the shared scratch, since the home directory is lost when the server stops.
default_base="$HOME/mesa"
[ -n "${JUPYTERHUB_USER:-}" ] && [ -w /user/scratch14_outreach ] && default_base="/user/scratch14_outreach/$JUPYTERHUB_USER/mesa"
base=${2:-$default_base}
row=$(echo "$TABLE" | awk -F'|' -v v="$want" '$1 == v')
[ -n "$row" ] || { echo "Unknown version '$1'. Known:"; "$0" --list; exit 1; }
IFS='|' read -r ver rec zip md5 folder sver srec sfile smd5 sfolder <<< "$row"

MESA_DIR="$base/$folder"
SDK_PARENT="$base/sdk-$sver"
MESASDK_ROOT="$SDK_PARENT/$sfolder"
mkdir -p "$base" "$base/downloads"
echo "MESA $ver with SDK $sver into $base"

fetch() {   # record file md5 -> verified file in $base/downloads
  local f="$base/downloads/$2"
  if [ -f "$f" ] && [ "$(md5sum "$f" | cut -d' ' -f1)" = "$3" ]; then echo "  have $2"; return; fi
  echo "  downloading $2"
  curl -fL --retry 3 -o "$f.part" "https://zenodo.org/records/$1/files/$2?download=1"
  [ "$(md5sum "$f.part" | cut -d' ' -f1)" = "$3" ] || { echo "md5 mismatch for $2"; rm -f "$f.part"; exit 1; }
  mv "$f.part" "$f"
}

if [ ! -f "$MESASDK_ROOT/bin/mesasdk_init.sh" ]; then
  fetch "$srec" "$sfile" "$smd5"
  mkdir -p "$SDK_PARENT"
  echo "  unpacking SDK $sver"
  tar -xzf "$base/downloads/$sfile" -C "$SDK_PARENT"
fi
if [ ! -d "$MESA_DIR" ]; then
  fetch "$rec" "$zip" "$md5"
  echo "  unpacking MESA $ver"
  unzip -q "$base/downloads/$zip" -d "$base"
fi

export MESASDK_ROOT MESA_DIR
export HDF5_USE_FILE_LOCKING=FALSE   # network file systems
set +u; . "$MESASDK_ROOT/bin/mesasdk_init.sh"; set -u   # first use also writes $MESASDK_ROOT/etc

# Known fixes for old versions on newer compilers (MESA FAQ); harmless where the pattern is absent.
case "$ver" in
  12115|10398)
    grep -rl "stop'fixup'" "$MESA_DIR/star/private" 2>/dev/null | xargs -r sed -i "s/stop'fixup'/stop 'fixup'/g" ;;
esac
if [ "$ver" = 10398 ] && [ -f "$MESA_DIR/star/private/history.f90" ]; then
  sed -i 's/if (open_close_log) close(io)/if (open_close_log .and. write_flag) close(io)/' "$MESA_DIR/star/private/history.f90"
fi

cd "$MESA_DIR"
echo "  compiling (log: $MESA_DIR/install.log; this takes a while)"
if ! ./install > install.log 2>&1; then
  echo "MESA $ver did not compile. The end of $MESA_DIR/install.log:"; tail -20 install.log; exit 1
fi
rm -f "$base/downloads/$zip" "$base/downloads/$sfile"   # 1-3 GB not needed once unpacked and compiled
echo
echo "MESA $ver is installed. To use it in a shell (or add these to ~/.bashrc):"
echo "  export MESASDK_ROOT=$MESASDK_ROOT"
echo "  source \$MESASDK_ROOT/bin/mesasdk_init.sh"
echo "  export MESA_DIR=$MESA_DIR"
