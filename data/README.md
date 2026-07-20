# Data folder

This folder is kept for local input files only. In the public R script, input paths are written as `/path/to/<file name>` placeholders so that private lab or cluster paths are not exposed.

Expected filenames:

- `inter_30.hic`
- `merged_loops.bedpe`
- `10000_blocks.bedpe`
- `CTCF.bed`
- `RAD21.liver.bed`
- `CtrlVeh_H3K27ac_log2.bw`
- `CtrlVeh_BRD4_log2.bw`
- `genes.gtf`

Large genomics files are ignored by `.gitignore` and should usually not be committed to a public GitHub repository. If the data are public, link to the data source in the main README. If the data are controlled, unpublished, or patient-derived, keep them outside the repository. Do not commit absolute paths that identify your private system.
