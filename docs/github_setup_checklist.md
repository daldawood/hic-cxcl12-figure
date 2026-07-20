# GitHub setup checklist for this project

## Before you make the repository public

- Replace `Your Name` in `LICENSE`.
- Replace author and repository fields in `CITATION.cff`.
- Remove private cluster paths, project names, accessions, unpublished sample names, and protected health information.
- Do not commit raw or large genomics files unless you have permission to share them publicly.
- Add a small test dataset only if it is public and redistributable.
- Add a figure preview to the README if the image can be shared.

## Suggested first commit

```bash
git init
git add README.md CITATION.cff LICENSE .gitignore data/README.md results/.gitkeep docs scripts
git commit -m "Add CXCL12 regional Hi-C figure workflow"
```

## Suggested GitHub repository settings

- Repository name: `hic-cxcl12-figure`
- Description: `R workflow for plotting a CXCL12-centered regional Hi-C figure with regulatory tracks`
- Topics: `bioinformatics`, `hic`, `genomics`, `r`, `ggplot2`, `cxcl12`, `chromatin`, `3d-genome`
- Visibility: Public only after checking data/privacy requirements.


## Privacy check before publishing

- Replace private paths with placeholders such as `/path/to/inter_30.hic`.
- Do not upload protected, unpublished, patient-derived, or very large genomics files.
- Search the repository for lab names, usernames, cluster names, and institution-specific folders before making it public.
