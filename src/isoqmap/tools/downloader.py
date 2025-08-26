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
    "gencode_38_v42": {
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
    "gencode_38": {
        "X_matrix": {
            "url": "https://github.com/ZhengCQ/IsoQMap/releases/download/v1.0.0/gencode_38.v39.X_matrix.RData.gz",
            "sha256": "445bce7cb49f11c08b505dfa00d2b6d9666142160e7f37c050476882ee19692c",
            "filename": "X_matrix.RData.gz"
        },
        "transcript": {
            "url": "https://github.com/ZhengCQ/IsoQMap/releases/download/v1.0.0/gencode_38.v39.transcript.fa.gz",
            "sha256": "7eb1097943b2169ff1f808350d96ae6f32593fd7fdd9e29f63e12ed7e5c7de81",
            "filename": "transcript.fa.gz"  
        },
        "geneinfo":{
            "url": "https://github.com/ZhengCQ/IsoQMap/releases/download/v1.0.0/gencode_38.v39.transcript_gene_info.tsv.gz",
            "sha256": "ce44c902780d03030e86c61af21a0fe7cb2a470181e8d67e3487d940b6ac667d",
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
            "url": "https://github.com/ZhengCQ/IsoQMap/releases/download/v1.0.0/refseq_38.110.transcript.fa.gz",
            "sha256": "61539d23a315690c0d5d609aef9c956c596abfd4e7d0206c6dabc75b56ceceb7",
            "filename": "transcript.fa.gz"
            
        },
         "geneinfo":{
            "url": "https://github.com/ZhengCQ/IsoQMap/releases/download/v1.0.0/refseq_38.110.transcript_gene_info.tsv.gz",
            "sha256": "c4cd130026914d72c866b9f1e11773ab3946f60336ae7efaf57fb9d8112e90fc",
            "filename": "transcript_gene_info.tsv.gz"
        } 
    },
    "pig_110":{
        "X_matrix": {
            "url": "https://github.com/ZhengCQ/IsoQMap/releases/download/v1.0.0/pig_110.X_matrix.RData.gz",
            "sha256": "900cd4a7e037e3ac11eb9b0d0c08f7b3fea488321a16b7d000d8312d647e5795",
            "filename": "X_matrix.RData.gz"
        },
        "transcript": {
            "url": "https://github.com/ZhengCQ/IsoQMap/releases/download/v1.0.0/pig_110.transcript.fa.gz",
            "sha256": "09379a4f747525eea821a1f56e79a6dacfe4a4a2f3f0c9d43e3fa1c8a37ed53d",
            "filename": "transcript.fa.gz"
        }, 
         "geneinfo":{
            "url": "https://github.com/ZhengCQ/IsoQMap/releases/download/v1.0.0/pig_110.transcript_gene_info.tsv.gz",
            "sha256": "70c520618afb11137e7ab072eb282e49d28e6e37ea1c4dd98696a98ac4fcb9ba",
            "filename": "transcript_gene_info.tsv.gz"
            }     
    }
}

RESOURCE_ROOT = Path(__file__).resolve().parent.parent / "resources" / "ref"

def decompress_zip(file_path):
    output_dir = Path(file_path).parent
    print(f"Decompressing {file_path} -> {output_dir}")
    with zipfile.ZipFile(file_path, 'r') as zip_ref:
        zip_ref.extractall(output_dir)
    print(f"✔ Decompressed to directory: {output_dir}")

# -------------------------------
# 计算文件 SHA256
# -------------------------------
def sha256sum(filepath):
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()

# -------------------------------
# 检查是否需要下载
# -------------------------------
def need_download(dest, expected_sha256, expected_size=None):
    if not dest.exists():
        return True
    if expected_size and dest.stat().st_size != expected_size:
        print(f"⚠ File size mismatch: {dest.stat().st_size} != {expected_size}")
        return True
    if sha256sum(dest) != expected_sha256:
        print(f"⚠ SHA256 mismatch for {dest}")
        return True
    return False

# -------------------------------
# 下载文件，支持断点续传
# -------------------------------
def download_file_with_resume(url, dest_path, retries=10, delay=5):
    part_path = dest_path.with_suffix(dest_path.suffix + ".part")

    for attempt in range(1, retries + 1):
        try:
            resume_byte_pos = os.path.getsize(part_path) if part_path.exists() else 0
            headers = {"Range": f"bytes={resume_byte_pos}-"} if resume_byte_pos > 0 else {}

            print(f"\nAttempt {attempt}/{retries}: Downloading {url} (resume from {resume_byte_pos} bytes)")

            with requests.get(url, headers=headers, stream=True, timeout=30) as response:
                if response.status_code in (200, 206):
                    total = int(response.headers.get("content-length", 0)) + (resume_byte_pos if response.status_code == 206 else 0)
                    mode = "ab" if response.status_code == 206 else "wb"
                    if response.status_code == 206:
                        print("✔ Server supports resume.")
                    else:
                        if resume_byte_pos > 0:
                            print("⚠ Server does not support resume, restarting from beginning...")
                            resume_byte_pos = 0

                    downloaded = resume_byte_pos
                    with open(part_path, mode) as f:
                        for chunk in response.iter_content(chunk_size=8192):
                            if chunk:
                                f.write(chunk)
                                downloaded += len(chunk)
                                if total:
                                    done = int(50 * downloaded / total)
                                    print(f"\r[{'█'*done}{'.'*(50-done)}] {downloaded/total:.2%}", end="")
                    print(f"\n✔ Download finished: {part_path}")

                    # 下载完成后重命名为正式文件
                    part_path.rename(dest_path)
                    return True
                else:
                    raise Exception(f"Unexpected status code: {response.status_code}")

        except KeyboardInterrupt:
            print(f"\n⏸ Download interrupted by user. Partial file saved as {part_path}")
            return False
        except Exception as e:
            print(f"\n✘ Download failed: {e}")
            if attempt < retries:
                print(f"Retrying after {delay} seconds...")
                time.sleep(delay)
            else:
                print("✘ Exceeded maximum retries. Download aborted.")
                return False

# -------------------------------
# 解压 .gz 文件
# -------------------------------
def decompress_gz(gz_path):
    dest_path = gz_path.with_suffix('')
    print(f"Decompressing {gz_path} -> {dest_path}")
    with gzip.open(gz_path, 'rb') as f_in, open(dest_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    print(f"✔ Decompressed: {dest_path}")
    return dest_path

# -------------------------------
# 下载参考文件
# -------------------------------
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

        if need_download(dest, meta["sha256"], meta.get("size")):
            success = download_file_with_resume(meta["url"], dest)
            if not success:
                raise RuntimeError(f"Failed to download: {filename}")

        # 下载完成后再次校验 SHA256
        if sha256sum(dest) != meta["sha256"]:
            print(f"✘ Hash mismatch after download. Deleting file.")
            dest.unlink()
            raise ValueError(f"Hash mismatch for {filename} after download.")

        print(f"✔ Downloaded and verified: {dest}")

        # 自动解压 .gz 文件（例如 X_matrix 文件）
        if name == "X_matrix" and dest.suffix == ".gz":
            decompress_gz(dest)

def download_osca():
    dest_dir = str(Path(__file__).resolve().parent.parent / "resources")
    download_file_with_resume('https://yanglab.westlake.edu.cn/software/osca/download/osca-0.46.1-linux-x86_64.zip',
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
