"""
Lightweight logging facade for BASALT-Air.

The legacy pipeline mixes ``print()`` with manual ``open('Basalt_log.txt', 'a')``
writes scattered across ~30 modules. This module wraps stdlib ``logging`` so
new code can call a single ``log(...)`` (or ``log.info``, ``log.warning`` etc.)
to emit timestamped messages to **both** stdout and the workdir log file
without further wiring. Existing prints/writes keep working unchanged.

Typical usage from the CLI / module entry points::

    from basalt.logger import setup_logger, get_logger
    setup_logger()              # call once after chdir to workdir
    log = get_logger()
    log.info("Pipeline started")
    log.warning("Something looks off")

If ``setup_logger`` is never called, ``get_logger()`` returns a logger that
is configured lazily on first use with sensible defaults pointing at
``Basalt_log.txt`` in the current working directory.
"""

from __future__ import annotations

import logging
import os
import sys

LOGGER_NAME = "basalt"
DEFAULT_LOGFILE = "Basalt_log.txt"

_DEFAULT_FORMAT = "[%(asctime)s] [%(levelname)s] %(message)s"
_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"

_initialised = False


def setup_logger(logfile=None, level=logging.INFO):
    """Configure the BASALT-Air logger to write to stdout and ``logfile``.

    Idempotent: subsequent calls re-target the file handler if ``logfile``
    differs (useful when the CLI chdirs into ``--workdir`` after parsing).
    """
    global _initialised
    logger = logging.getLogger(LOGGER_NAME)
    logger.setLevel(level)
    logger.propagate = False

    target = os.path.abspath(logfile) if logfile else os.path.abspath(DEFAULT_LOGFILE)

    have_stream = False
    for h in list(logger.handlers):
        if isinstance(h, logging.StreamHandler) and not isinstance(h, logging.FileHandler):
            have_stream = True
        elif isinstance(h, logging.FileHandler):
            if os.path.abspath(getattr(h, 'baseFilename', '')) == target:
                # Already pointing at the right file; keep it.
                continue
            logger.removeHandler(h)
            h.close()

    formatter = logging.Formatter(_DEFAULT_FORMAT, _DATE_FORMAT)

    if not have_stream:
        sh = logging.StreamHandler(stream=sys.stdout)
        sh.setFormatter(formatter)
        logger.addHandler(sh)

    file_attached = any(
        isinstance(h, logging.FileHandler)
        and os.path.abspath(getattr(h, 'baseFilename', '')) == target
        for h in logger.handlers
    )
    if not file_attached:
        try:
            fh = logging.FileHandler(target, mode='a', encoding='utf-8')
            fh.setFormatter(formatter)
            logger.addHandler(fh)
        except OSError as e:
            # Don't crash the pipeline just because the log file is unwritable.
            sys.stderr.write(
                "[BASALT-Air] WARNING: could not open log file {!r}: {}\n".format(target, e)
            )

    _initialised = True
    return logger


def get_logger():
    """Return the singleton BASALT-Air logger, configuring it lazily if needed."""
    if not _initialised:
        setup_logger()
    return logging.getLogger(LOGGER_NAME)


def format_elapsed(seconds):
    """Render a duration in seconds as a human-friendly H:MM:SS string."""
    seconds = max(0, int(round(float(seconds))))
    hours, rem = divmod(seconds, 3600)
    minutes, sec = divmod(rem, 60)
    return "{:d}:{:02d}:{:02d}".format(hours, minutes, sec)
