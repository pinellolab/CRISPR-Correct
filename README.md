# CRISPR-Correct
**CRISPR-Correct** is a Python package for performing guide RNA mapping given the raw FASTQs and guide library dataframe. Specifically, <ins>CRISPR-Correct is able to handle imperfect mapping </ins> (due to self-editing using SpRY-based base-editors or sequencing errors) by mapping observed protospacer sequences to the guide RNA with the closest hamming distance; CRISPR-Correct is unable to handle in-dels in the protospacer.

For several large samples, CRISPR-Correct can also be ran on Broad Institute's Terra Platform. The workflow file is located at: https://portal.firecloud.org/?return=terra#methods/pinellolab/CrisprSelfEditMappingOrchestratorWorkflowSampleEntity/2

