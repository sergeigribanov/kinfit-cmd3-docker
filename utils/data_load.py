import os
from tqdm.auto import tqdm
import shutil
import tarfile
import json
import requests
import urllib.parse as ul

def load_tar(url, prefix='.'):
    full_url = 'https://cloud-api.yandex.net:443/v1/disk/public/resources/download?public_key=' + url
    res = json.loads(requests.get(full_url).text)
    filename = ul.parse_qs(ul.urlparse(res['href']).query)['filename'][0]
    tar_path = os.path.join(prefix, filename)
    with requests.get(res['href'], stream=True) as r:
        total_length = int(r.headers.get("Content-Length"))
        with tqdm.wrapattr(r.raw, "read", total=total_length, desc="")as raw:
            with open(tar_path, 'wb')as output:
                shutil.copyfileobj(raw, output)

    with tarfile.open(tar_path, 'r') as fl:
        fl.extractall(prefix)

    os.remove(tar_path)
