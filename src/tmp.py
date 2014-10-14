# remove soft-clipping for STAR->RSEM

# awk 'BEGIN {OFS="\t"} {split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}' Aligned.out.sam > Aligned.noS.sam

import sys, re
out_sep='\t'
number_pat = re.compile('[0-9]*')
smdin_pat = re.compile('[SMDIN]')
numberS_pat = re.compile('[0-9]*S')
CIGAR_I = 5
SEQ_I = 9
QUAL_I = 10
for row in sys.stdin:
    fields = row.split(' ')
    # split on integers and on [SMDIN]
    c = re.findall(number_pat, fields[MAPQUAL_I])
    l = re.findall(smdin_pat, fields[MAPQUAL_I])
    # hocus pocus to remove soft clipping
    if c[2] == 'S':
        fields[SEQ_I]=fields[SEQ_I][l[1]+1:]
        fields[QUAL_I]=fields[QUAL_I][l[1]+1:]
    if c[-1] == 'S':
        l1=len(fields[SEQ_I]) - l[len(l)-1]
        fields[SEQ_I]=fields[SEQ_I][1:l1]
        fields[MAPQUAL_I]=re.sub(numberS_pat, '', fields[MAPQUAL_I])
    sys.stdout.write(out_sep.join(row))
