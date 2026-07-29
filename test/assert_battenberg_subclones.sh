#!/bin/bash
function md5_battenberg_subclones {
    ### get columns 1 to 12: chr startpos endpos BAF pval LogR ntot nMaj1_A nMin1_A frac1_A nMaj2_A nMin2_A
    cat "$1" | cut -f 1-12 | md5sum | cut -f 1 -d ' '
}

received=$(md5_battenberg_subclones "$1")
expected=$(md5_battenberg_subclones "$2")

if [ "$received" == "$expected" ]; then
    echo "Battenberg subclones files are equal"
    exit 0
else
    echo "Battenberg subclones files files are not equal" >&2
    exit 1
fi
