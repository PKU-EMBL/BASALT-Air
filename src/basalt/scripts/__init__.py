"""Bundled non-Python scripts (Perl helpers).

Use :func:`path` to resolve a script name to its on-disk location after
``pip install``::

    from basalt.scripts import path as _script_path
    subprocess.run(['perl', _script_path('calc.kmerfreq.pl'), ...])
"""

from importlib.resources import files


def path(name):
    """Return the absolute filesystem path of a bundled script.

    Parameters
    ----------
    name : str
        Filename of the bundled resource (e.g. ``'calc.kmerfreq.pl'``).
    """
    return str(files(__package__) / name)
