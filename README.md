# kentoabeywardane.github.io

Personal website — a single static page. No build step, no dependencies.

```
index.html              # all content
assets/css/style.css    # all styling (light + dark via prefers-color-scheme)
assets/img/*.webp       # project figures
.nojekyll               # serve files as-is on GitHub Pages
```

## Preview locally

```bash
python3 -m http.server 8000   # then open http://localhost:8000
```

## Deploy on GitHub Pages

1. Create a public repo named **`kentoabeywardane.github.io`**.
2. `git remote add origin git@github.com:kentoabeywardane/kentoabeywardane.github.io.git && git push -u origin main`
3. Settings → Pages → Source: *Deploy from a branch*, branch `main`, folder `/ (root)`.

The site goes live at `https://kentoabeywardane.github.io` within a minute or two.

## Editing

Everything is plain HTML in `index.html`, in document order: intro, featured research,
publications, other projects, experience. Colors, spacing, and type live in the
token block at the top of `style.css`.

To add a CV: drop the PDF at `assets/Kento_Abeywardane_CV.pdf` and uncomment the CV
link in the `.links` list in `index.html`.
