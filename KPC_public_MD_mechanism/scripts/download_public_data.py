"""Download and verify the public inputs used by the KPC analysis.

Only compact MD-derived feature/energy archives and six PDB coordinate files are
downloaded. The much larger QM/MM coordinate-path archives are not required.
"""

from __future__ import annotations

import argparse
import hashlib
import tarfile
import urllib.request
from pathlib import Path

# Checksums are the values published by Zenodo for record 7114981.
ZENODO_FILES = {
    "2.enes.tar.gz": "2bdb00cb6bd463713359a911dec28a11",
    "3.datasets.tar.gz": "76f646c858d2946135b00f152a9c129e",
}
PDB_IDS = ("5UL8", "7TB7", "7TBX", "7TC1", "4ZBE", "8AKL")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/raw/public_kpc"),
        help="Directory for archives, extracted arrays, and PDB files",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Download again even when a local file passes validation",
    )
    return parser.parse_args()


def md5sum(path: Path) -> str:
    """Return the published-style MD5 checksum for a downloaded archive."""
    digest = hashlib.md5()  # noqa: S324 -- integrity check, not cryptographic security
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, destination: Path) -> None:
    """Download to a temporary file, then atomically replace the destination."""
    temporary = destination.with_suffix(destination.suffix + ".part")
    request = urllib.request.Request(url, headers={"User-Agent": "KPC-public-data/1.0"})
    with urllib.request.urlopen(request, timeout=120) as response, temporary.open("wb") as out:
        while chunk := response.read(1024 * 1024):
            out.write(chunk)
    temporary.replace(destination)


def safe_extract(archive: Path, output_dir: Path) -> None:
    """Extract a tar archive only when every member remains below output_dir."""
    root = output_dir.resolve()
    with tarfile.open(archive, "r:gz") as handle:
        for member in handle.getmembers():
            if member.issym() or member.islnk():
                raise ValueError(f"Archive links are not allowed: {member.name}")
            target = (output_dir / member.name).resolve()
            if root not in target.parents and target != root:
                raise ValueError(f"Unsafe archive member: {member.name}")
        handle.extractall(output_dir)  # noqa: S202 -- paths were checked above


def fetch_zenodo(output_dir: Path, force: bool) -> None:
    """Download, checksum, and extract the two compact Zenodo archives."""
    for name, expected_md5 in ZENODO_FILES.items():
        destination = output_dir / name
        valid = destination.exists() and md5sum(destination) == expected_md5
        if force or not valid:
            url = f"https://zenodo.org/records/7114981/files/{name}?download=1"
            print(f"Downloading {url}")
            download(url, destination)
        observed = md5sum(destination)
        if observed != expected_md5:
            destination.unlink(missing_ok=True)
            raise ValueError(f"Checksum mismatch for {name}: {observed}")
        safe_extract(destination, output_dir)


def fetch_pdb_files(output_dir: Path, force: bool) -> None:
    """Download PDB coordinate files and reject obvious HTML/error responses."""
    for pdb_id in PDB_IDS:
        destination = output_dir / f"{pdb_id}.pdb"
        valid = destination.exists() and "\nATOM  " in destination.read_text(
            encoding="utf-8", errors="ignore"
        )
        if force or not valid:
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            print(f"Downloading {url}")
            download(url, destination)
        text = destination.read_text(encoding="utf-8", errors="ignore")
        if "\nATOM  " not in text:
            destination.unlink(missing_ok=True)
            raise ValueError(f"Downloaded file does not look like a PDB: {pdb_id}")


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    fetch_zenodo(args.output_dir, args.force)
    fetch_pdb_files(args.output_dir, args.force)
    print(f"Public KPC inputs are ready in {args.output_dir}")


if __name__ == "__main__":
    main()
