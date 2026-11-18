#!/usr/bin/env python3
"""Pack a ShEPhERD trajectory (ChimeraX .xyz/.bild frame dump) into a compact,
base64-encoded JS payload for the web viewer.

Usage:
  python3 tools/build_traj.py <frame_dir> <out.js> [--var NAME] [--stride N]

Input layout (as written by the ShEPhERD trajectory dumper):
  NNN_frame_x1.xyz   atoms           (element + xyz)
  NNN_frame_x2.bild  surface points  (unused; geometry duplicated by x3)
  NNN_frame_x3.bild  ESP surface     (.color encodes signed ESP: blue=+, red=-)
  NNN_frame_x4.bild  pharmacophores  (.sphere anchors, .arrow direction vectors)
"""
import argparse, base64, glob, json, os, re, sys
import numpy as np

POS_SCALE = 819.0  # int16 full scale -> +/- 40 A, ~0.0012 A resolution


def read_xyz(path):
    lines = [l for l in open(path).read().splitlines() if l.strip()]
    n = int(lines[0])
    rows = lines[1:]
    if len(rows) == n + 1:          # optional comment line
        rows = rows[1:]
    el, xyz = [], []
    for r in rows[:n]:
        t = r.split()
        el.append(t[0])
        xyz.append([float(v) for v in t[1:4]])
    return el, np.array(xyz, dtype=np.float64).reshape(-1, 3)


def read_bild(path):
    """-> (spheres, arrows); sphere=(color, xyz, r); arrow=(color, start, end)."""
    col = (0.0, 0.0, 0.0)
    spheres, arrows = [], []
    for line in open(path):
        t = line.split()
        if not t:
            continue
        if t[0] == '.color':
            col = tuple(float(v) for v in t[1:4])
        elif t[0] == '.sphere':
            spheres.append((col, [float(v) for v in t[1:4]], float(t[4])))
        elif t[0] == '.arrow':
            arrows.append((col, [float(v) for v in t[1:4]], [float(v) for v in t[4:7]]))
    return spheres, arrows


def b64(arr):
    return base64.b64encode(np.ascontiguousarray(arr).tobytes()).decode('ascii')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('frame_dir')
    ap.add_argument('out')
    ap.add_argument('--var', default='SHEPHERD2_TRAJ')
    ap.add_argument('--stride', type=int, default=3)
    a = ap.parse_args()

    idx = sorted(int(re.search(r'(\d+)_frame_x1\.xyz$', p).group(1))
                 for p in glob.glob(os.path.join(a.frame_dir, '*_frame_x1.xyz')))
    if not idx:
        sys.exit(f'no frames in {a.frame_dir}')
    # Always keep the final frame, whatever the stride.
    keep = sorted(set(idx[::a.stride]) | {idx[-1]})
    print(f'{len(idx)} frames -> {len(keep)} kept (stride {a.stride})')

    # ---- pass 1: read everything, collect the element + colour palettes -------
    frames = []
    elements, colors = [], []

    def cidx(c):
        key = tuple(round(v, 4) for v in c)
        if key not in colors:
            colors.append(key)
        return colors.index(key)

    for i in keep:
        s = f'{i:03d}'
        el, p1 = read_xyz(os.path.join(a.frame_dir, f'{s}_frame_x1.xyz'))
        for e in el:
            if e not in elements:
                elements.append(e)
        sph3, _ = read_bild(os.path.join(a.frame_dir, f'{s}_frame_x3.bild'))
        sph4, arr4 = read_bild(os.path.join(a.frame_dir, f'{s}_frame_x4.bild'))
        frames.append(dict(el=el, p1=p1, s3=sph3, s4=sph4, a4=arr4))

    # ---- centre on the final frame's molecule, so playback settles centred ----
    center = frames[-1]['p1'].mean(axis=0)

    def q(p):
        return np.round((np.asarray(p, dtype=np.float64) - center) * POS_SCALE).astype(np.int16)

    x1_counts, x1_elem, x1_pos = [], [], []
    x3_pos, x3_esp = [], []
    x4s_counts, x4s_col, x4s_pos, x4s_rad = [], [], [], []
    x4a_counts, x4a_col, x4a_pos, x4a_dir, x4a_len = [], [], [], [], []

    for f in frames:
        x1_counts.append(len(f['el']))
        x1_elem.extend(elements.index(e) for e in f['el'])
        x1_pos.append(q(f['p1']))

        # x3: .color is (max(-esp,0), 0, max(esp,0)) -- recover the signed field.
        # ChimeraX clamps channel values at 1.0, so clip at +/-1 in raw units to
        # keep the same saturation the reference rendering shows.
        pos3 = np.array([s[1] for s in f['s3']])
        esp = np.array([s[0][2] - s[0][0] for s in f['s3']])
        x3_pos.append(q(pos3))
        x3_esp.append(np.round(np.clip(esp, -1, 1) * 127).astype(np.int8))

        x4s_counts.append(len(f['s4']))
        for c, p, r in f['s4']:
            x4s_col.append(cidx(c))
            x4s_pos.append(q(p))
            x4s_rad.append(int(round(r * 50)))

        x4a_counts.append(len(f['a4']))
        for c, p0, p1_ in f['a4']:
            v = np.array(p1_) - np.array(p0)
            L = np.linalg.norm(v) or 1.0
            x4a_col.append(cidx(c))
            x4a_pos.append(q(p0))
            x4a_dir.append(np.round(v / L * 127).astype(np.int8))
            x4a_len.append(int(min(255, round(L * 20))))

    payload = {
        'nframes': len(frames),
        'srcFrames': len(idx),
        'frameIndex': keep,
        'posScale': POS_SCALE,
        'elements': elements,
        'colors': [list(c) for c in colors],
        'x1': {'counts': b64(np.array(x1_counts, np.uint16)),
               'elem': b64(np.array(x1_elem, np.uint8)),
               'pos': b64(np.concatenate(x1_pos))},
        'x3': {'n': len(frames[0]['s3']),
               'pos': b64(np.concatenate(x3_pos)),
               'esp': b64(np.concatenate(x3_esp))},
        'x4s': {'counts': b64(np.array(x4s_counts, np.uint16)),
                'col': b64(np.array(x4s_col, np.uint8)),
                'rad': b64(np.array(x4s_rad, np.uint8)),
                'pos': b64(np.array(x4s_pos, np.int16).reshape(-1, 3))},
        'x4a': {'counts': b64(np.array(x4a_counts, np.uint16)),
                'col': b64(np.array(x4a_col, np.uint8)),
                'len': b64(np.array(x4a_len, np.uint8)),
                'pos': b64(np.array(x4a_pos, np.int16).reshape(-1, 3)),
                'dir': b64(np.array(x4a_dir, np.int8).reshape(-1, 3))},
    }

    with open(a.out, 'w') as fh:
        fh.write(f'/* generated by tools/build_traj.py from {os.path.basename(a.frame_dir)} */\n')
        fh.write(f'window.{a.var} = ')
        json.dump(payload, fh, separators=(',', ':'))
        fh.write(';\n')
    print(f'wrote {a.out}  ({os.path.getsize(a.out)/1024:.0f} KB)')
    print('element palette:', elements)
    print('colour palette :', colors)


if __name__ == '__main__':
    main()
