"""Scutoid census over all 31 PlantSeg ovule volumes: download files (4 at a time, since OSF is slow per
connection), run scutoid.py on each as it arrives, delete it. Finished volumes (data/scutoid/<name>.json) are
skipped, so the run can be restarted."""
import concurrent.futures as cf
import pathlib
import subprocess

from run_all import NODES, PY, RAW, download, listing

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "data" / "scutoid"

if __name__ == "__main__":
    todo = []
    for node, (kind, _, vox) in NODES.items():
        if kind == "ovules":
            todo += [(name, url, vox) for name, url in listing(node)
                     if not (OUT / (pathlib.Path(name).stem + ".json")).exists()]
    print(len(todo), "volumes to do", flush=True)

    def fetch(item):
        name, url, vox = item
        download(url, RAW / name)
        return item

    with cf.ThreadPoolExecutor(4) as ex:
        for name, url, vox in (f.result() for f in cf.as_completed([ex.submit(fetch, it) for it in todo])):
            dst = RAW / name
            r = subprocess.run([PY, str(HERE / "scutoid.py"), str(dst), *map(str, vox)], capture_output=True, text=True)
            print((r.stdout or r.stderr).strip().splitlines()[-1], flush=True)
            dst.unlink()
