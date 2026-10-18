#!/bin/sh
# Download OpenAI's proof sources (lean/OAI/Geometry/PeriodicTiling/*.lean) at the pinned commit into ./OAI/.
set -e
REV=adc7f1241b42e322a6451854ab7e4b4c146bf78a
cd "$(dirname "$0")"
mkdir -p OAI/Geometry/PeriodicTiling
curl -sS "https://api.github.com/repos/openai/math/git/trees/$REV?recursive=1" |
  python3 -c "import sys,json; [print(t['path']) for t in json.load(sys.stdin)['tree'] if t['path'].startswith('lean/OAI/Geometry/PeriodicTiling/') and t['path'].endswith('.lean')]" |
  while read -r p; do curl -sS "https://raw.githubusercontent.com/openai/math/$REV/$p" -o "${p#lean/}"; done
ls OAI/Geometry/PeriodicTiling/*.lean | wc -l
