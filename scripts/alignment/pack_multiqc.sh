#!/bin/bash
source activate snakemake-rnaseq

multiqc fastp_reports/ aligned_star/ --title "Preeclampsia placenta RNA-seq report" -b ""

tar -czf multiqc_report.tar.gz Preeclampsia_placenta_RNA-seq_report_multiqc_report.html \ 
             Preeclampsia_placenta_RNA-seq_report_multiqc_report/
