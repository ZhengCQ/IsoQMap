import hashlib
import urllib.request
from pathlib import Path
import requests
import time
import os
import gzip
import shutil
import zipfile


# 定义支持的数据版本及其对应文件和哈希值
REFERENCE_DATA = {
    "gencode_38": {
        "X_matrix": {
            "url": "https://github.com/ZhengCQ/IsoQMap/releases/download/v1.0.0/gencode_38.v41.X_matrix.RData.gz",
            "sha256": "f088e4f29e9d582fca4e6e4b46a7e08a8358d89a3661c910bbe73c44a80e52d0",
            "filename": "X_matrix.RData.gz"
        },
        "transcript": {
            "url": "https://github.com/ZhengCQ/IsoQMap/releases/download/v1.0.0/gencode_38.v41.transcript.fa.gz",
            "sha256": "172d04be1deaf2fd203c2d9063b2e09b33e3036dd2f169d57d996a6e8448fe94",
            "filename": "transcript.fa.gz"  
        },
        "geneinfo":{
            "url": "https://github.com/ZhengCQ/IsoQMap/releases/download/v1.0.0/gencode_38.v41.transcript_gene_info.tsv.gz",
            "sha256": "f93ed5707479af4072d26a324b9193a348d403878d93823c9cbf933a59d6261c",
            "filename": "transcript_gene_info.tsv.gz"
            } 
    },
    "refseq_38": {
        "X_matrix": {
            "url": "https://github.com/ZhengCQ/IsoQMap/releases/download/v1.0.0/refseq_38.110.X_matrix.RData",
            "sha256": "9c758d2177065e0d8ae4fc8b5d6bcb3d45e7fe8f9a0151669a1eee230f2992d1",
            "filename": "X_matrix.RData.gz"

        },
        "transcript": {
            "url": "https://github.com/ZhengCQ/IsoQMap/releases/download/v1.0.0/refseq_38.transcript.fa.gz",
            "sha256": "..."
        },
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
}

RESOURCE_ROOT = Path(__file__).resolve().parent.parent / "resources" / "ref"

def sha256sum(file_path):
    h = hashlib.sha256()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()

def decompress_gz(file_path):
    output_path = file_path.with_suffix('')
    print(f"Decompressing {file_path} -> {output_path}")
    with gzip.open(file_path, 'rb') as f_in, open(output_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    print(f"✔ Decompressed: {output_path}")
    
def decompress_zip(file_path):
    output_dir = Path(file_path).parent
    print(f"Decompressing {file_path} -> {output_dir}")
    with zipfile.ZipFile(file_path, 'r') as zip_ref:
        zip_ref.extractall(output_dir)
    print(f"✔ Decompressed to directory: {output_dir}")

def download_file_with_retry(url, dest_path, retries=6, delay=3):
    for attempt in range(1, retries + 1):
        try:
            # 1. 检查本地是否已有部分内容
            resume_byte_pos = os.path.getsize(dest_path) if os.path.exists(dest_path) else 0

            headers = {"Range": f"bytes={resume_byte_pos}-"} if resume_byte_pos > 0 else {}
            print(f"Attempt {attempt} to download {url} (resuming from {resume_byte_pos} bytes)")

            with requests.get(url, headers=headers, stream=True, timeout=10) as response:
                if response.status_code not in [200, 206]:
                    raise Exception(f"Unexpected status code: {response.status_code}")
                
                total = int(response.headers.get('content-length', 0)) + resume_byte_pos
                downloaded = resume_byte_pos

                mode = 'ab' if resume_byte_pos > 0 else 'wb'
                with open(dest_path, mode) as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
                            downloaded += len(chunk)
                            done = int(50 * downloaded / total) if total else 0
                            print(f"\r[{'█' * done}{'.' * (50 - done)}] {downloaded / total:.2%}", end='')

            print(f"\n✔ Download succeeded: {dest_path}")
            return True

        except Exception as e:
            print(f"\n✘ Download failed: {e}")
            if attempt < retries:
                print(f"Retrying after {delay} seconds...")
                time.sleep(delay)
            else:
                print("✘ Exceeded retry limit.")
                return False


def download_reference(version='gencode_38', files_requested=['all']):
    if version not in REFERENCE_DATA:
        raise ValueError(f"Unsupported reference version: {version}")

    version_dir = RESOURCE_ROOT / version
    version_dir.mkdir(parents=True, exist_ok=True)

    for name, meta in REFERENCE_DATA[version].items():
        if 'all' not in files_requested and name not in files_requested:
            continue

        filename = meta["filename"]
        dest = version_dir / filename

        if dest.exists():
            print(f"{dest} already exists. Verifying hash...")
            if sha256sum(dest) == meta["sha256"]:
                print(f"✔ Hash OK for {filename}, skipping download.")
            else:
                print(f"✘ Hash mismatch for {filename}, re-downloading...")
                dest.unlink()
        if not dest.exists():
            success = download_file_with_retry(meta["url"], dest)
            if not success:
                raise RuntimeError(f"Failed to download: {filename}")

        print(f"Verifying downloaded file {filename}...")
        if sha256sum(dest) != meta["sha256"]:
            print(f"✘ Hash mismatch after download. Deleting file.")
            dest.unlink()
            raise ValueError(f"Hash mismatch for {filename} after download.")

        print(f"✔ Downloaded and verified: {dest}")

        # 自动解压X_matrix文件（如果是.gz结尾）
        if name == "X_matrix" and dest.suffix == ".gz":
            decompress_gz(dest)

def download_osca():
    dest_dir = str(Path(__file__).resolve().parent.parent / "resources")
    download_file_with_retry('https://yanglab.westlake.edu.cn/software/osca/download/osca-0.46.1-linux-x86_64.zip',
                            dest_dir + '/' + 'osca-0.46.1-linux-x86_64.zip')
    decompress_zip(dest_dir + '/' + 'osca-0.46.1-linux-x86_64.zip')
    os.system(f'chmod 755 {dest_dir}/osca-0.46.1-linux-x86_64/osca && ln -fs {dest_dir}/osca-0.46.1-linux-x86_64/osca {dest_dir}/osca')
    return f'{dest_dir}/osca' 

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Download reference resources")
    parser.add_argument("version", choices=REFERENCE_DATA.keys(), help="Reference version to download")
    parser.add_argument("--files", default="all", help="Comma-separated file types to download (default: all)")
    args = parser.parse_args()
    files_requested = args.files.split(',') if args.files else ['all']

    download_reference(args.version, files_requested)
