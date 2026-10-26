"""01_download.py —— 下载 GEO 数据（Python 版）

为什么需要它：本机 shell 的 Windows 凭证存储被沙箱限制，
PowerShell / curl / R 的 download.file 都会在 TLS 握手阶段失败；
Python 用自带 OpenSSL，且本机出口有代理，故需跳过证书校验
（下载的是 NCBI 公开数据，可接受）。

用法： python R/01_download.py
"""
import os
import ssl
import sys
import urllib.request

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(PROJ, "data")
os.makedirs(DATA, exist_ok=True)

B421 = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217421/suppl/"
B423 = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217423/suppl/"
B490 = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253490/suppl/"

MUST = [
    (B421, "GSE217421_Conv-RNAseq-Read-Counts.tsv.gz", 10),        # MB
    (B421, "GSE217421_Conv-RNAseq-Experiment-Design.tsv.gz", 0.001),
    (B421, "GSE217421_Program-Source-Codes.tar.gz", 0.01),
    (B421, "GSE217421_README.pdf.gz", 0.01),
]

CTX = ssl._create_unverified_context()


def fetch(base, name, min_mb):
    dest = os.path.join(DATA, name)
    if os.path.exists(dest) and os.path.getsize(dest) / 1e6 >= min_mb:
        print(f"[skip] {name}  ({os.path.getsize(dest)/1e6:.1f} MB)")
        return
    url = base + name
    print(f"[get ] {name}")
    with urllib.request.urlopen(url, context=CTX, timeout=300) as r, open(dest, "wb") as f:
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            f.write(chunk)
    print(f"       -> {os.path.getsize(dest)/1e6:.2f} MB")


if __name__ == "__main__":
    for base, name, mb in MUST:
        fetch(base, name, mb)
    print("[done] 数据目录:", DATA)
