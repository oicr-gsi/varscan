#!/bin/bash
cd $1


echo ".vcf files:"
for v in *.vcf;do gunzip -c $v | grep -v GATKCommandLine | md5sum;done | sort -V

echo ".indel files:"
find . -name *.indel | xargs md5sum | sort -V

echo ".snp files:"
find . -name *.snp | xargs md5sum | sort -V

echo ".copynumber files:"
find . -name *.copynumber | xargs md5sum | sort -V

echo ".copynumber.filtered files:"
find . -name *.copynumber.filtered | xargs md5sum | sort -V
