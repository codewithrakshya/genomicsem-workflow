# Publishing this repository

This repository contains the reusable pipeline and fictional examples. Review
the staged file list before committing changes:

```bash
git status --short
git diff --cached --stat
```

The repository uses the MIT license in `LICENSE`. Dependencies retain their own
licenses and citation requirements. The GitHub repository is:
https://github.com/codewithrakshya/genomicsem-workflow

For future updates:

```bash
git add PATHS_YOU_REVIEWED
git diff --cached
git commit -m "Describe the update"
git push origin main
```

Keep raw data, reference downloads, results, run manifests, and local R libraries
out of commits. `.gitignore` protects the standard locations, but does not prevent
force-adds or data placed elsewhere. Use config/local/ for private configuration.
Avoid replacing the tracked fictional examples with identifiable study details.
