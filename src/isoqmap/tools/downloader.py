import hashlib
import urllib.request
from pathlib import Path
import requests
import time
import os

# 定义支持的数据版本及其对应文件和哈希值
REFERENCE_DATA = {
    "gencode_38": {
        "X_matrix": {
            "url": "https://github.com/ZhengCQ/IsoQMap/releases/download/v1.0.0/gencode_38.X_matrix.RData.gz",
            "sha256": "f088e4f29e9d582fca4e6e4b46a7e08a8358d89a3661c910bbe73c44a80e52d0"  
        },
        "transcript": {
            "url": "https://github.com/ZhengCQ/IsoQMap/releases/download/v1.0.0/gencode_38.transcript.fa.gz",
            "sha256": "172d04be1deaf2fd203c2d9063b2e09b33e3036dd2f169d57d996a6e8448fe94"  
        }, 
    },
    "refseq_38": {
        "X_matrix": {
            "url": "https://github.com/ZhengCQ/IsoQMap/releases/download/v1.0.0/refseq_38.110.X_matrix.RData",
            "sha256": "9c758d2177065e0d8ae4fc8b5d6bcb3d45e7fe8f9a0151669a1eee230f2992d1"
        },
        "transcript": {
            "url": "https://github.com/ZhengCQ/IsoQMap/releases/download/v1.0.0/refseq_38.transcript.fa.gz",
            "sha256": "..."
        },
    "pig_110":{
        "X_matrix": {
            "url": "https://github.com/ZhengCQ/IsoQMap/releases/download/v1.0.0/pig_110.X_matrix.RData.gz",
            "sha256": "900cd4a7e037e3ac11eb9b0d0c08f7b3fea488321a16b7d000d8312d647e5795"  
        },
        "transcript": {
            "url": "https://github.com/ZhengCQ/IsoQMap/releases/download/v1.0.0/pig_110.transcript.fa.gz",
            "sha256": "09379a4f747525eea821a1f56e79a6dacfe4a4a2f3f0c9d43e3fa1c8a37ed53d"  
        }, 
        
    }
    },
    # 可添加更多版本
}



# 设置本地保存路径
RESOURCE_ROOT = Path(__file__).resolve().parent.parent / "resources" / "ref"

def sha256sum(file_path):
    h = hashlib.sha256()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def download_file_with_retry(url, dest_path, retries=3, delay=2):
    for attempt in range(1, retries + 1):
        try:
            print(f"Attempt {attempt} to download {url}")
            response = requests.get(url, stream=True, timeout=10)
            response.raise_for_status()

            with open(dest_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            print(f"Download succeeded. Save in the {dest_path}")
            return True
        except Exception as e:
            print(f"Download failed: {e}")
            if os.path.exists(dest_path):
                os.remove(dest_path)
            if attempt < retries:
                time.sleep(delay)
            else:
                print("Exceeded retry limit.")
                return False


def download_reference(version='gencode_38', files_requested=['all']):
    if version not in REFERENCE_DATA:
        raise ValueError(f"Unsupported reference version: {version}")

    version_dir = RESOURCE_ROOT / version
    version_dir.mkdir(parents=True, exist_ok=True)

    for name, meta in REFERENCE_DATA[version].items():
        if 'all' not in files_requested and name not in files_requested:
            continue

        filename = meta["url"].split("/")[-1].replace(version + '.', '')
        dest = version_dir / filename

        if dest.exists():
            print(f"{dest} already exists. Verifying hash...")
            if sha256sum(dest) == meta["sha256"]:
                print(f"✔ Hash OK for {filename}, skipping download.")
                continue
            else:
                print(f"✘ Hash mismatch for {filename}, re-downloading...")

        download_file_with_retry(meta["url"], dest)
        print(f"Verifying downloaded file {filename}...")
        if sha256sum(dest) != meta["sha256"]:
            raise ValueError(f"Hash mismatch for {filename} after download.")

        print(f"✔ Downloaded and verified: {dest}")




if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Download reference resources")
    parser.add_argument("version", choices=REFERENCE_DATA.keys(), help="Reference version to download")
    args = parser.parse_args()
    files_requested = args.files.split(',') if args.files else ['all']

    download_reference(args.version, files_requested)

