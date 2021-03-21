#!/bin/bash
source activate snakemake-rnaseq

cd fastp_reports

find . -name '*_fastp_report.json' -exec bash -c 'cp $0 ${0/_report/}' {} \;

cd ..

multiqc fastp_reports/ aligned_star/ --title "Preeclampsia placenta RNA-seq report" -b "" -v

tar -czf multiqc_report.tar.gz Preeclampsia-placenta-RNA-seq-report_multiqc_report.html Preeclampsia-placenta-RNA-seq-report_multiqc_report_data/

