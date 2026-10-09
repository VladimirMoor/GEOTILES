#!/bin/zsh
# Перепись полидрафтеров: числа Хееша (heesch-sat, Kaplan), быстрый прогон до 5 слоёв.
# Использование: census.sh n1 n2 ...   (результаты: data/drafterN_out.txt, сводка data/drafterN_report.txt)
SRC=/Users/vladimirmuravev/Projects2026/GEOTILES/tools/heesch-sat/src
D=/Users/vladimirmuravev/Projects2026/GEOTILES/research/p3/data
for n in "$@"; do
  [ -f $D/drafter$n.txt ] || $SRC/gen -drafter -size $n -free -o $D/drafter$n.txt
  ( start=$(date +%s); $SRC/sat -isohedral -maxlevel 5 $D/drafter$n.txt -o $D/drafter${n}_out.txt;
    $SRC/report $D/drafter${n}_out.txt > $D/drafter${n}_report.txt; echo "time $(( $(date +%s) - start ))s" >> $D/drafter${n}_report.txt ) &
done
wait
