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
2. `git remote add origin https://github.com/kentoabeywardane/kentoabeywardane.github.io.git && git push -u origin master`
3. Settings → Pages → Source: *Deploy from a branch*, branch `master`, folder `/ (root)`.

The site goes live at `https://kentoabeywardane.github.io` within a minute or two.

Editing notes, the interactive viewers, and how to rebuild their data are in
[`assets/NOTES.md`](assets/NOTES.md).
