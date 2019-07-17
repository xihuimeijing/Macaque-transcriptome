#!/bin/sh
less -S statis.L1Single.iden75.Topres |sed -n '6,16p'|cut -f 1,2 > box.finalfigure75 
Rscript /mnt/share/lisx/scripts/barplot.R -i=box.finalfigure75 -x="CDS coverage " -c=grey -y="percent" -o="barplot-final.identity75.pdf" -m="CDS protein coverage (identity: 75%)"

less -S statis.all.iden75.Topres |sed -n '6,16p'|cut -f 1,2 > box.finalfigureall75
Rscript barplot.R -i=box.finalfigureall75 -x="CDS coverage " -c=grey -y="percent" -o="barplot-final-alllist.identity75.pdf" -m="CDS protein coverage (identity: 75%)"




