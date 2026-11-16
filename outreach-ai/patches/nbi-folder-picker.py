#!/usr/bin/env python3
"""Notebook Intelligence 6.3.0: the start-folder picker becomes a text box (card 1113-1).

NBI's launcher tiles for Claude Code and OpenCode ask for a start folder with JupyterLab's FileDialog. Pressing Enter after
naming a new folder there closes the dialog before the rename: "File browser is disposed", the folder stays "Untitled Folder".
This build step replaces the one async function Ln (the picker) in the minified chunk 128 with a text-input dialog that
creates the folder through the contents API. It returns the same value as before: a path relative to the server root, no
leading slash, or undefined when cancelled.

The server sends chunks with "Cache-Control: public, max-age=31536000, immutable", so a browser that has the old chunk never
asks again. The patched chunk therefore gets a new file name (128.<20 hex of its sha256>.js), the remoteEntry that names it
gets a new name too (remoteEntry.<20 hex>.js), and package.json's jupyterlab._build.load points at the new remoteEntry.
The old two files are deleted.

Guarded: every pinned sha256 (chunk, remoteEntry, package.json) and every expected occurrence count is checked on every file
before anything is written. On any mismatch it prints why and exits 1, so a different NBI version, or a second run, fails
the image build instead of getting a patch written for another one. Standard library only. New files get the owner, group
and mode of the files they replace; package.json is rewritten in place.

Usage: nbi-folder-picker.py [labextension dir ...]   (no argument: the two installed copies; arguments are for testing)
"""
import hashlib
import json
import os
import sys

# The extension is installed twice: where JupyterLab serves it, and inside the Python package it was copied from.
ROOTS = [
    "/opt/conda/share/jupyter/labextensions/@plmbr/notebook-intelligence",
    "/opt/conda/lib/python3.13/site-packages/notebook_intelligence/labextension",
]
OLD_CHUNK_HASH = "060938e2eb7426187b40"
OLD_CHUNK = "static/128.%s.js" % OLD_CHUNK_HASH
OLD_REMOTE_HASH = "a38b1fbd478d0e6e12b9"
OLD_REMOTE = "static/remoteEntry.%s.js" % OLD_REMOTE_HASH
PIN = {
    OLD_CHUNK: "060938e2eb7426187b40424bfff996f49c7169c797048df6543e4d311bf45e7a",
    OLD_REMOTE: "a38b1fbd478d0e6e12b90f29dfac3ebf857cd10dfebee23c103a167b144f487b",
    "package.json": "167a06855e84c797d142aa3501e3a1ac27ff0356436b54aeb86026e201987585",
}
MARKER = "/*outreach-ai patch: nbi-folder-picker, card 1113-1*/"

ORIGINAL = (
    'async function Ln(e,t,n){const o=await at.FileDialog.getExistingDirectory({manager:e,title:t,label:t,defaultPath:n});'
    'if(!o.button.accept)return;const a=o.value;return a&&a.length>0?a[0].path:null!=n?n:""}'
)

# e = IDocumentManager (e.services is the ServiceManager), t = dialog title, n = current file-browser folder (no leading
# slash, "" at the root). `a` is the module's @jupyterlab/apputils. A leading ~ means /home/user. Contains no single quote,
# so it fits the r'''...''' strings below.
REPLACEMENT = (
    r'''async function Ln(e,t,n){''' + MARKER + r'''try{'''
    r'''const o=null!=n&&""!==String(n).replace(/^\/+/,"")?"/"+String(n).replace(/^\/+/,""):"/home/user";'''
    # The help line under the title is Falk's wording (card 1113-1, "K").
    r'''const i=await a.InputDialog.getText({title:t,label:"Type a new path to make a folder, or open the folder in the file browser first; this box starts there.",text:o,okLabel:"OK"});'''
    r'''if(!i.button.accept)return;'''
    r'''const s=String(null!=i.value?i.value:"").trim(),l=(""===s?o:s),r=("~"===l||l.startsWith("~/")?"/home/user"+l.slice(1):l).split("/").filter(e=>""!==e&&"."!==e);'''
    r'''if(r.includes("..")){await a.showErrorMessage(t,"The folder path must not contain .. segments.");return}'''
    r'''const c=(e.services||e.serviceManager).contents;'''
    r'''let u="";'''
    r'''for(const e of r){u=""===u?e:u+"/"+e;let n=null;'''
    r'''try{n=await c.get(u,{content:!1})}catch(e){if(!e||!e.response||404!==e.response.status)throw e}'''
    r'''if(n&&"directory"!==n.type){await a.showErrorMessage(t,"/"+u+" exists and is a file, not a folder.");return}'''
    r'''if(!n)await c.save(u,{type:"directory"})}'''
    r'''return u}catch(o){try{await a.showErrorMessage(t,String(o&&o.message||o))}catch(e){}}}'''
)


def fail(msg):
    print("nbi-folder-picker: FAILED: " + msg)
    sys.exit(1)


def sha(b):
    return hashlib.sha256(b).hexdigest()


def read(path):
    if not os.path.isfile(path):
        fail("%s is missing" % path)
    with open(path, "rb") as f:
        return f.read()


def check(root):
    """Check every guard on one extension folder; return what is needed to patch it."""
    files = {}
    for rel, pinned in PIN.items():
        p = os.path.join(root, rel)
        d = read(p)
        if sha(d) != pinned:
            fail("%s has sha256 %s, expected %s (a different Notebook Intelligence build, or already patched)"
                 % (p, sha(d), pinned))
        files[rel] = d
    chunk = files[OLD_CHUNK].decode("utf-8")
    if chunk.count(ORIGINAL) != 1:
        fail("%s: the original Ln definition occurs %d times, expected exactly 1" % (root, chunk.count(ORIGINAL)))
    remote = files[OLD_REMOTE].decode("utf-8")
    # The chunk hash occurs twice in remoteEntry, both times as the entry for chunk 128: the file-name map and the ?v= map.
    if remote.count(OLD_CHUNK_HASH) != 2 or remote.count('128:"%s"' % OLD_CHUNK_HASH) != 2:
        fail("%s: the chunk hash occurs %d times in remoteEntry, expected 2 (as 128:\"<hash>\" in the name map and the ?v= map)"
             % (root, remote.count(OLD_CHUNK_HASH)))
    pkg_text = files["package.json"].decode("utf-8")
    old_load = '"load": "%s"' % OLD_REMOTE
    if pkg_text.count(old_load) != 1:
        fail("%s: package.json names the old remoteEntry %d times, expected 1" % (root, pkg_text.count(old_load)))
    pkg = json.loads(pkg_text)
    if pkg.get("jupyterlab", {}).get("_build", {}).get("load") != OLD_REMOTE:
        fail("%s: package.json jupyterlab._build.load is not %s" % (root, OLD_REMOTE))

    new_chunk = chunk.replace(ORIGINAL, REPLACEMENT).encode("utf-8")
    if new_chunk.count(MARKER.encode()) != 1:
        fail("%s: marker count after patching is not 1" % root)
    new_hash = sha(new_chunk)[:20]
    new_remote = remote.replace(OLD_CHUNK_HASH, new_hash).encode("utf-8")
    new_remote_hash = sha(new_remote)[:20]
    new_chunk_rel = "static/128.%s.js" % new_hash
    new_remote_rel = "static/remoteEntry.%s.js" % new_remote_hash
    new_pkg_text = pkg_text.replace(old_load, '"load": "%s"' % new_remote_rel)
    new_pkg = json.loads(new_pkg_text)
    # Only that one key may differ.
    expect = json.loads(pkg_text)
    expect["jupyterlab"]["_build"]["load"] = new_remote_rel
    if new_pkg != expect:
        fail("%s: package.json changed in more than the load key" % root)
    for rel in (new_chunk_rel, new_remote_rel):
        if os.path.exists(os.path.join(root, rel)):
            fail("%s already exists" % os.path.join(root, rel))
    return {
        "root": root, "new_chunk": (new_chunk_rel, new_chunk), "new_remote": (new_remote_rel, new_remote),
        "new_pkg": new_pkg_text.encode("utf-8"), "old_sha": sha(files[OLD_CHUNK]),
    }


def put_new(root, rel, data, like_rel):
    """Write a new file with the owner, group and mode of the file it replaces."""
    p = os.path.join(root, rel)
    st = os.stat(os.path.join(root, like_rel))
    with open(p, "xb") as f:
        f.write(data)
    os.chown(p, st.st_uid, st.st_gid)
    os.chmod(p, st.st_mode & 0o7777)
    with open(p, "rb") as f:
        if f.read() != data:
            fail("%s does not read back as written" % p)


def apply(plan):
    root = plan["root"]
    new_chunk_rel, new_chunk = plan["new_chunk"]
    new_remote_rel, new_remote = plan["new_remote"]
    put_new(root, new_chunk_rel, new_chunk, OLD_CHUNK)
    put_new(root, new_remote_rel, new_remote, OLD_REMOTE)
    pp = os.path.join(root, "package.json")
    with open(pp, "r+b") as f:  # in place: same inode, so owner, group, mode and ACL stay
        f.seek(0)
        f.write(plan["new_pkg"])
        f.truncate()
    os.remove(os.path.join(root, OLD_CHUNK))
    os.remove(os.path.join(root, OLD_REMOTE))
    print("nbi-folder-picker: patched %s -> %s (%s)" % (plan["old_sha"][:12], sha(new_chunk)[:12], root))
    print("nbi-folder-picker: chunk %s, remoteEntry %s, package.json load updated" % (new_chunk_rel, new_remote_rel))


def main():
    roots = sys.argv[1:] or ROOTS
    plans = [check(r) for r in roots]  # every guard on every file before anything is written
    if len({p["new_chunk"][0] for p in plans}) != 1 or len({p["new_remote"][0] for p in plans}) != 1:
        fail("the extension copies would get different file names")
    for plan in plans:
        apply(plan)


if __name__ == "__main__":
    main()
