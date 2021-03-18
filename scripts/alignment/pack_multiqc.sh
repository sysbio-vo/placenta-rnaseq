#!/bin/bash
source activate snakemake-rnaseq

multiqc fastp_reports/ aligned_star/ --title "Preeclampsia placenta RNA-seq report" -b "" -v

tar -czf multiqc_report.tar.gz Preeclampsia-placenta-RNA-seq-report_multiqc_report.html Preeclampsia-placenta-RNA-seq-report_multiqc_report_data/

