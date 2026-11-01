# Uploading this package to GitHub

1. Choose and add the software license you want before making the repository public.
2. Create an empty GitHub repository, for example `Complex-Weinberg-3HDM-BFB`.
3. From the root of this package run:

```bash
git init
git add .
git commit -m "Initial public release"
git branch -M main
git remote add origin git@github.com:YOUR-ACCOUNT/Complex-Weinberg-3HDM-BFB.git
git push -u origin main
```

If you use GitHub CLI, the create-and-push step can instead be:

```bash
gh repo create Complex-Weinberg-3HDM-BFB --public --source=. --remote=origin --push
```

After the paper has an arXiv identifier or DOI, update `CITATION.cff` and the citation section of `README.md`. Recommended repository topics include `3hdm`, `higgs`, `bounded-from-below`, `particle-physics`, `symbolic-computation`, and `sympy`.
