#!/usr/bin/env python

"""
Utility script for downloading and unpacking BASALT model weights.

This script can be invoked directly or from within the BASALT pipeline
to ensure that the required deep learning models are present locally.
"""

import os
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
    local_path = f"{local_dir}/BASALT.zip"

    if os.path.exists(local_path):
        print(f"File already exists at {local_path}.")
        return local_path

    os.makedirs(local_dir, exist_ok=True)
    print(f"File will be saved in {local_path}.")

    response = requests.get(url, stream=True)
    response.raise_for_status()
    total_size = int(response.headers.get('content-length', 0))
    block_size = 1024
    progress_bar = tqdm(total=total_size, unit='iB', unit_scale=True)

    with open(local_path, 'wb') as f:
        for data in response.iter_content(block_size):
            progress_bar.update(len(data))
            f.write(data)

    progress_bar.close()
    print("File downloaded successfully.")
    return local_path


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

    exit_code = os.system(f"unzip -o {zip_path} -d {destination_folder}")

    if exit_code == 0:
        print("Successfully unzipped BASALT model archive.")
    else:
        print("Unzip failed: please check your unzip configuration.")


if __name__ == '__main__':
    main()
