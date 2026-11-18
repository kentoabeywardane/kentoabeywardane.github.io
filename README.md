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

Plain HTML in `index.html`, in document order: profile sidebar, About Me, News,
Featured Work, Publications, Other Projects, Background. Colors, spacing, and type
live in the token block at the top of `style.css`.

- **Photo:** save one as `assets/img/portrait.jpg` and swap the `<div class="avatar monogram">`
  for the commented-out `<img class="avatar" ...>` right above it.
- **CV:** drop the PDF at `assets/Kento_Abeywardane_CV.pdf` and uncomment the CV link
  in `.profile-links`.
- **News:** newest first; keep it to about five items.

## Interactive trajectory viewer

`shepherd2-trajectory.html` plays back precomputed ShEPhERD reverse-diffusion
trajectories in the browser (three.js from a CDN; no build step).

### Adding or changing a trajectory

1. Pack a ChimeraX frame dump (the `NNN_frame_x{1,3,4}.{xyz,bild}` files that
   `shepherd_traj.py` loads) into a JS payload:

   ```sh
   python3 tools/build_traj.py <frame_dir> assets/data/traj-<name>.js \
       --var TRAJ_<NAME> --stride 3
   ```

   `--stride 3` keeps every third frame (~190 KB for 360 steps). Raise it for a
   smaller file; the last frame is always kept.

2. Add an entry to the `SCENES` array in the `<script type="module">` block at
   the bottom of `shepherd2-trajectory.html`:

   ```js
   { id: 'name', label: 'Tab label',
     src: 'assets/data/traj-name.js', global: 'TRAJ_NAME',
     note: 'One line shown under the tabs.' },
   ```

   **The array order is the tab order.** Payloads load lazily, so only the
   selected trajectory is fetched.

### Local preview

The viewer is an ES module, which browsers refuse to load over `file://`. Serve
the directory instead:

```sh
python3 -m http.server 8787
# then open http://127.0.0.1:8787/shepherd2-trajectory.html
```

`index.html` still opens fine directly from disk. Opening
`shepherd2-trajectory.html` from disk shows a short note saying to serve it
instead, rather than sitting on "Loading trajectory…".

This only affects local preview — on GitHub Pages everything is served over
https and the viewer works normally.
