#!/bin/sh
/mnt/share/lisx/tools/DIP/dpi1/identify_tss_peaks.sh -g /share/data/chr.size/rheMac8.size -i 'all-ctss.bed.gz/*.bed.gz' -o DIP-out 2> DIP.log 
less DIP-out/outPooled/tc.spi_merged.ctssMaxCounts11_ctssMaxTpm1.bed |awk '{if ($6=="+"){print $1"\t"$2"\t"$2+1}else{ print $1"\t"$3-1"\t"$3 }}' > peak.position
