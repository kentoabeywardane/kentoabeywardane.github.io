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

## Interactive viewers

`shepherd2-trajectory.html` carries two WebGL viewers (three.js from a CDN, no
build step):

| | |
|---|---|
| `assets/js/evolution-viewer.js` | the seven-stage fragment-merging story |
| `assets/js/traj-viewer.js` | one trajectory on its own, further down the page |

Shared pieces: `mol-core.js` (decoding + instanced molecule drawing),
`ribbon.js` (protein cartoon), `viewpoint.js` (choosing an opening camera).

### The fragment-merging story

Everything lives in the 8CNX crystal frame, translated by one shared offset so
the world origin is the binding site. Nothing is re-centred per structure —
that is what makes the stage-to-stage overlays meaningful.

```sh
python3 tools/build_evolution.py \
    <.../shepherd2_composition_fragmerge_8-10_edm_trajectory> --stride 3
```

Writes three payloads into `assets/data/`:

- `evo-base.js` (~23 KB) — the protein CA trace plus the five static
  structures (fragments 8/10/12, the round-1 merge, the round-2 product).
  Loaded as a blocking script, since the first frame needs it.
- `evo-traj-1.js`, `evo-traj-2.js` (~65 KB each) — the two generation
  trajectories, fetched the first time their stage is opened. Atoms only:
  the x3 and x4 channels denoise alongside them, but a half-noisy surface
  reads as nothing and hides the molecule underneath it, so they are not
  packed. The profiles that are shown are the fragments', which are a fixed
  target and live in `evo-base.js`.

The protein comes from `assets/data/8cnx_apo.pdb`. Only CA positions, the
carbonyl direction that orients the ribbon, and a secondary-structure flag
leave Python; the mesh is built in the browser.

**Stages** are the `STAGES` array at the top of `evolution-viewer.js` — a
`lead` list (solid, may show its interaction profile), a `refs` list (ghosted
context), whether the profile starts visible, and the caption. `CARBON` just
above it maps a structure to its carbon colour. Adding or reordering a stage
means editing that array and nothing else.

`#stage-4` in the URL opens the page straight on that stage, with no
transition in from stage 1.

### A single trajectory

```sh
python3 tools/build_traj.py <frame_dir> assets/data/traj-<name>.js \
    --var TRAJ_<NAME> --stride 3
```

Unlike `build_evolution.py`, this re-centres on the final frame, which is what
you want for a trajectory shown on its own. `--stride 3` keeps every third
frame (~190 KB for 360 steps); the last frame is always kept.

### shepherd-score alignment

`shepherd-score-alignment.html` replays one single-start shepherd-score
alignment in three representations (volume without H, surface, pharmacophores),
with tabs to switch between them. Viewer: `assets/js/align-viewer.js`. Payload:
`assets/data/align-82-99.js` (~44 KB), built from the trajectory dump that
`gen_alignment_trajs.py` writes (one `mol<i>__mol<j>/` folder):

```sh
# needs rdkit + shepherd_score (the shepherdscore conda env)
python tools/build_alignment.py <.../mol82__mol99> \
    <.../moses_test_scaffold_molblock_charges.pkl> assets/data/align-82-99.js
```

The fit moves rigidly, so the payload stores it once in its start pose, plus
one quaternion and translation per optimizer step (recovered by Kabsch from the
dumped frames). Pharmacophores are recomputed from the molblocks; surface
points are read from the BILD files, since sampling them is random.

### Local preview

The viewers are ES modules, which browsers refuse to load over `file://`.
Serve the directory instead:

```sh
python3 -m http.server 8787
# then open http://127.0.0.1:8787/shepherd2-trajectory.html
```

`index.html` still opens fine directly from disk. Opening the viewer page from
disk shows a short note saying to serve it instead, rather than hanging on
"Loading…". This only affects local preview — on GitHub Pages everything is
served over https.
