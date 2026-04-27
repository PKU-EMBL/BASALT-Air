#!/usr/bin/env python

"""
Lightweight subprocess helper for BASALT.

Use ``run_cmd`` instead of ``os.system`` for any new or critical external
tool invocation: it surfaces non-zero exit codes, captures stderr for
diagnostics, and supports timeouts so a stuck subprocess cannot hang the
pipeline indefinitely.
"""

import shlex
import subprocess
import sys


class CommandError(RuntimeError):
    """Raised when an external command exits non-zero or times out."""


def run_cmd(cmd, timeout=None, check=True, capture=False, cwd=None, env=None):
    """
    Run an external command with proper error reporting.

    Parameters
    ----------
    cmd : str or list of str
        Command to run. Strings are split with ``shlex.split`` and executed
        without a shell, which avoids accidental glob/redirect injection.
        Pass a list directly for arguments that contain spaces.
    timeout : float or None
        Kill the process and raise ``CommandError`` if it runs longer than
        this many seconds. ``None`` disables the timeout (not recommended
        for long-running external tools).
    check : bool
        If True, raise ``CommandError`` on non-zero exit. If False, return
        the ``CompletedProcess`` regardless of exit code.
    capture : bool
        If True, capture stdout and stderr into the returned object.
        Otherwise they stream to the caller's terminal (default, matches
        ``os.system`` behaviour).
    cwd : str or None
        Working directory for the child process.
    env : dict or None
        Environment for the child process (defaults to the parent's).

    Returns
    -------
    subprocess.CompletedProcess

    Raises
    ------
    CommandError
        If the command exits non-zero (and ``check=True``) or exceeds
        ``timeout``. The message includes the command, exit code, and the
        last chunk of stderr when available.
    """
    if isinstance(cmd, str):
        argv = shlex.split(cmd)
        display = cmd
    else:
        argv = list(cmd)
        display = ' '.join(shlex.quote(str(a)) for a in argv)

    stdout = subprocess.PIPE if capture else None
    stderr = subprocess.PIPE if capture else None

    try:
        result = subprocess.run(
            argv,
            stdout=stdout, stderr=stderr,
            timeout=timeout, cwd=cwd, env=env,
            text=True,
        )
    except subprocess.TimeoutExpired as e:
        raise CommandError(
            'command timed out after {}s: {}'.format(timeout, display)
        ) from e
    except FileNotFoundError as e:
        raise CommandError(
            'executable not found for command: {} ({})'.format(display, e)
        ) from e

    if check and result.returncode != 0:
        tail = (result.stderr or '').strip().splitlines()[-20:]
        tail_text = '\n  '.join(tail) if tail else '<no stderr captured>'
        raise CommandError(
            'command failed (exit {}): {}\n  stderr tail:\n  {}'.format(
                result.returncode, display, tail_text
            )
        )
    return result


def run_shell(cmd, timeout=None, check=True, cwd=None, env=None):
    """
    Drop-in replacement for ``os.system`` that uses a shell.

    Prefer ``run_cmd`` when you do not actually need shell features
    (pipes, redirects, globs, ``&&``). This helper exists for migrating
    legacy call sites that rely on shell semantics without rewriting
    them. Return code behaviour mirrors ``run_cmd``.
    """
    try:
        result = subprocess.run(
            cmd,
            shell=True,
            timeout=timeout, cwd=cwd, env=env,
            text=True,
        )
    except subprocess.TimeoutExpired as e:
        raise CommandError(
            'shell command timed out after {}s: {}'.format(timeout, cmd)
        ) from e

    if check and result.returncode != 0:
        raise CommandError(
            'shell command failed (exit {}): {}'.format(result.returncode, cmd)
        )
    return result


if __name__ == '__main__':
    # Smoke test: should print "hello" and exit 0
    run_cmd(['echo', 'hello'])
    # Should raise CommandError
    try:
        run_cmd(['false'])
    except CommandError as e:
        sys.stderr.write('expected failure caught: {}\n'.format(e))
