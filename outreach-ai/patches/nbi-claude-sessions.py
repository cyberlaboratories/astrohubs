#!/usr/bin/env python3
"""Notebook Intelligence 6.3.0: the Claude Code launcher tile lists earlier sessions whatever the chat panel's mode (card 1113-6).

The launcher tile "Claude Code" opens a dialog that first asks GET /notebook-intelligence/claude-sessions for earlier sessions.
ClaudeSessionsListHandler.get answers 404 "Claude Code mode is not enabled" when ai_service_manager.is_claude_code_mode is false.
That flag is the CHAT panel's setting (claude_settings.enabled), and on this hub the chat panel deliberately uses open-weight
models, so the flag is off and the dialog showed a red error even though "+ New Session" starts the claude CLI in a terminal.
Picking a session in the tile opens a terminal with `cd <cwd> && claude --resume <id>`; it makes no server call.

This build step deletes that one 5-line guard from the list handler and puts a comment line in its place. The resume handler
(POST /claude-sessions/resume, it reconnects the chat client) has an identical guard and is left alone: the snippet that is
matched starts at "class ClaudeSessionsListHandler", so the other guard cannot be matched.

Guarded: the full sha256 of extension.py is pinned and the snippet must occur exactly once. On any mismatch it prints why and
exits 1, so a different NBI version, or a second run, fails the image build. Standard library only. The file is rewritten in
place (same inode: owner, group, mode and ACL stay).

Usage: nbi-claude-sessions.py [extension.py]   (no argument: the installed copy; an argument is for testing)
"""
import hashlib
import os
import py_compile
import sys
import tempfile

TARGET = "/opt/conda/lib/python3.13/site-packages/notebook_intelligence/extension.py"
PIN = "371272d5bb39cd50ad6745da6a3256184d6ada3642e48a58a4d76fc1d1ee4a6f"
MARKER = "outreach-ai patch: nbi-claude-sessions, card 1113-6"

HEAD = '''class ClaudeSessionsListHandler(APIHandler):
    """Lists Claude Code sessions for both pickers.

    Both pickers (chat sidebar and launcher tile) consume this endpoint
    via the same ``list_all_sessions`` backend so they cannot disagree on
    previews. ``?scope=cwd`` filters the response down to sessions whose
    transcript lives under the current Jupyter cwd (the chat sidebar
    case); ``?scope=all`` (the default) returns every project's sessions
    (the launcher tile case).

    ``current_cwd`` is the realpath-resolved Jupyter root the chat picker
    pairs with the session id to render ``cd ... && claude --resume <id>``.

    Guarded by ``is_claude_code_mode`` because both consumers are
    Claude-Code-specific surfaces — the launcher tile literally launches
    the ``claude`` CLI in a terminal.
    """

    @tornado.web.authenticated
    def get(self):
'''
GUARD = '''        if not ai_service_manager.is_claude_code_mode:
            self.set_status(404)
            trans = i18n.request_bundle(self)
            self.finish(json.dumps({"error": trans.__("Claude Code mode is not enabled")}))
            return
'''
ORIGINAL = HEAD + GUARD
REPLACEMENT = HEAD + "        # " + MARKER + ": list sessions whatever the chat panel's mode\n"


def fail(msg):
    print("nbi-claude-sessions: FAILED: " + msg)
    sys.exit(1)


def sha(b):
    return hashlib.sha256(b).hexdigest()


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else TARGET
    if not os.path.isfile(path):
        fail("%s is missing" % path)
    with open(path, "rb") as f:
        old = f.read()
    if sha(old) != PIN:
        fail("%s has sha256 %s, expected %s (a different Notebook Intelligence build, or already patched)"
             % (path, sha(old), PIN))
    text = old.decode("utf-8")
    n = text.count(ORIGINAL)
    if n != 1:
        fail("the list handler's guard snippet occurs %d times, expected exactly 1" % n)
    guard_before = text.count("Claude Code mode is not enabled")
    new_text = text.replace(ORIGINAL, REPLACEMENT)
    new = new_text.encode("utf-8")
    if new_text.count(MARKER) != 1:
        fail("marker count after patching is not 1")
    if new_text.count("Claude Code mode is not enabled") != guard_before - 1:
        fail("the guard message count did not drop by exactly 1 (was %d)" % guard_before)
    # Syntax check without leaving a .pyc behind.
    with tempfile.TemporaryDirectory() as d:
        src = os.path.join(d, "extension.py")
        with open(src, "wb") as f:
            f.write(new)
        try:
            py_compile.compile(src, cfile=os.path.join(d, "extension.pyc"), doraise=True)
        except py_compile.PyCompileError as e:
            fail("the patched file does not compile: %s" % e)
    with open(path, "r+b") as f:  # in place: same inode, so owner, group, mode and ACL stay
        f.seek(0)
        f.write(new)
        f.truncate()
    with open(path, "rb") as f:
        if f.read() != new:
            fail("%s does not read back as written" % path)
    print("nbi-claude-sessions: patched %s -> %s" % (sha(old)[:12], sha(new)[:12]))


if __name__ == "__main__":
    main()
