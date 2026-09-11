#!/usr/bin/env python

"""
Utility script for downloading and unpacking BASALT model weights.

This script can be invoked directly or from within the BASALT-Air pipeline
to ensure that the required deep learning models are present locally.
"""

import os
import stat
import tempfile
import zipfile
from pathlib import Path

import requests
from tqdm import tqdm


def download_model(local_dir=None):
    """
    Download the BASALT model archive into the given directory.

    If ``local_dir`` is not provided, the directory pointed to by the
    BASALT_WEIGHT environment variable is used.
    """
    if local_dir is None:
        local_dir = os.environ.get("BASALT_WEIGHT")
        if not local_dir:
            raise EnvironmentError(
                "BASALT_WEIGHT environment variable is not set. "
                "Please configure a directory for BASALT model weights."
            )

    # Public URL hosting BASALT model weights
    url = "https://figshare.com/ndownloader/files/41093033"
    local_path = os.path.join(local_dir, "BASALT.zip")

    if os.path.exists(local_path):
        print(f"File already exists at {local_path}.")
        return local_path

    os.makedirs(local_dir, exist_ok=True)
    print(f"File will be saved in {local_path}.")

    # Download to a temporary file in the destination directory so an
    # interrupted transfer cannot leave a corrupt archive that looks complete.
    fd, temp_path = tempfile.mkstemp(prefix='.BASALT.', suffix='.part', dir=local_dir)
    os.close(fd)
    try:
        response = requests.get(url, stream=True, timeout=(10, 120))
        response.raise_for_status()
        total_size = int(response.headers.get('content-length', 0))
        block_size = 1024 * 1024
        with tqdm(total=total_size, unit='iB', unit_scale=True) as progress_bar:
            with open(temp_path, 'wb') as f:
                for data in response.iter_content(block_size):
                    if data:
                        progress_bar.update(len(data))
                        f.write(data)
        os.replace(temp_path, local_path)
    except Exception:
        try:
            os.unlink(temp_path)
        except OSError:
            pass
        raise
    print("File downloaded successfully.")
    return local_path


def _safe_extract(zip_path, destination):
    """Extract a model archive without allowing path traversal."""
    destination = Path(destination).resolve()
    destination.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path) as archive:
        for member in archive.infolist():
            mode = (member.external_attr >> 16) & 0xFFFF
            if stat.S_ISLNK(mode):
                raise ValueError(
                    'refusing to extract symbolic link from model archive: {}'
                    .format(member.filename)
                )
            target = (destination / member.filename).resolve()
            if target != destination and destination not in target.parents:
                raise ValueError(
                    'refusing to extract archive member outside destination: {}'
                    .format(member.filename)
                )
        archive.extractall(destination)


def main():
    """
    Download and unzip BASALT model weights into BASALT_WEIGHT.
    """
    zip_path = download_model()
    destination_folder = os.environ.get("BASALT_WEIGHT")
    if not destination_folder:
        raise EnvironmentError(
            "BASALT_WEIGHT environment variable is not set. "
            "Please configure a directory for BASALT model weights."
        )

    _safe_extract(zip_path, destination_folder)
    print("Successfully unzipped BASALT model archive.")


if __name__ == '__main__':
    main()
