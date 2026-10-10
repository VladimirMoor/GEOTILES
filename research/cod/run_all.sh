#!/bin/zsh
# Последовательный прогон по группам (один процесс, nice). Лог: data/run_all.log
for sg in "I 41 2 2" "P 61 2 2" "P 65 2 2" "I -4 2 d" "P 43 3 2" "P 41 3 2" "I 41 c d" "I 21 3" "P 21 3" "P 41 21 2" "P 43 21 2" "I 41" "P 31 2 1"; do
  nice -n 10 ../../.venv/bin/python cod_dv.py "$sg" 2>&1 | grep -v "^  [0-9]"
done
echo ALL DONE
