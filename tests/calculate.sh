#!/bin/bash
cd $1

echo ".vcf files:"
for v in *.vcf;do grep -v ^# $v | cut -f 1-6 | md5sum;done | sort -V

echo ".indel files:"
find . -name *.indel | cut -f 1-6 | xargs md5sum | sort -V

echo ".snp files:"
find . -name *.snp | cut -f 1-6 | xargs md5sum | sort -V

echo ".copynumber files:"
find . -name *.copynumber | cut -f 1-6 | xargs md5sum | sort -V

echo ".copynumber.filtered files:"
find . -name *.copynumber.filtered | cut -f 1-6 | xargs md5sum | sort -V
