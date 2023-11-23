# Somatic variant calling

This directory contains Nextflow code which implements the somatic variant calling pipeline used in the paper. 
Please note that this is not intended to run on new data and is instead intended to provide information on software, parameters, and workflow used. The 
software and parameters can be used independently without needing to use the Nextflow implementation. The pipeline can be reproduced by creating a Nextflow 
config for your given local computing environment and by editing the somatic_variant_pipeline.nf paths to reflect your local system directory.

The Singularity container containing the software versions used in the pipeline are hosted on [Zenodo](https://zenodo.org/records/10180501) due to 
Github file size limit constraints. BED files are included here and on Zenodo.

# Acknowledgements

We would like to thank David Marron, Alan Hoyle, and the UNC Lineberger Bioinformatics Core for establishing this pipeline and
assisting us in adapting it to our specific needs with BAM slices for this project.
