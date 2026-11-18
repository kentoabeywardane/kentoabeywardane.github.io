#!/usr/bin/env python3
"""Pack the bioisosteric fragment-merging story into payloads for the web viewer.

The story is told in one shared coordinate frame -- the 8CNX crystal frame --
because the whole point is that the fragments, the generated molecules and the
protein all sit on top of each other.  Everything is translated by a single
offset (the pocket centre) so the viewer can keep the world origin at the
binding site, and nothing is re-centred per structure the way a standalone
trajectory would be.

Writes:
  assets/data/evo-base.js   protein ribbon + the five static structures
  assets/data/evo-traj-1.js round-1 generation trajectory (frag 8 + 10)
  assets/data/evo-traj-2.js round-2 generation trajectory (+ frag 12)

Usage:
  python3 tools/build_evolution.py <merge_dir> [--stride 3] [--out assets/data]

<merge_dir> is the shepherd2_composition_fragmerge_8-10_edm_trajectory folder.
"""
import argparse, base64, glob, json, os, re, sys
import numpy as np

POS_SCALE = 819.0   # int16 -> +/- 40 A for ligand-scale coordinates
PRO_SCALE = 256.0   # int16 -> +/- 128 A; the protein is far bigger than the pocket


# --------------------------------------------------------------------------
# readers
# --------------------------------------------------------------------------
def read_xyz(path):
    lines = [l for l in open(path).read().splitlines() if l.strip()]
    n = int(lines[0])
    rows = lines[1:]
    if len(rows) == n + 1:                    # optional comment line
        rows = rows[1:]
    el, xyz = [], []
    for r in rows[:n]:
        t = r.split()
        el.append(t[0])
        xyz.append([float(v) for v in t[1:4]])
    return el, np.array(xyz, dtype=np.float64).reshape(-1, 3), None


def read_sdf(path):
    """-> (elements, coords, bonds).  V2000 only, which is what RDKit writes here."""
    L = open(path).read().splitlines()
    na, nb = int(L[3][0:3]), int(L[3][3:6])
    el, xyz = [], []
    for i in range(na):
        r = L[4 + i]
        xyz.append([float(r[0:10]), float(r[10:20]), float(r[20:30])])
        el.append(r[31:34].strip())
    bonds = []
    for i in range(nb):
        r = L[4 + na + i]
        bonds.append((int(r[0:3]) - 1, int(r[3:6]) - 1))
    return el, np.array(xyz, dtype=np.float64), bonds


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


# --------------------------------------------------------------------------
# a shared palette across every structure in the story
#
# Each structure used to carry its own element and colour tables.  Here they
# all share one, so a pharmacophore type is the same swatch in every stage and
# the viewer can hold a single legend.
# --------------------------------------------------------------------------
class Palette:
    def __init__(self):
        self.elements, self.colors = [], []

    def elem(self, e):
        if e not in self.elements:
            self.elements.append(e)
        return self.elements.index(e)

    def color(self, c):
        key = tuple(round(v, 4) for v in c)
        if key not in self.colors:
            self.colors.append(key)
        return self.colors.index(key)


def pack_structure(el, xyz, bonds, s3, s4, a4, center, pal):
    """One frame, in the layout the viewer reads for both statics and frames.

    Pass s3 as None to pack the atoms alone: the profile channels are then
    left out of the payload entirely and the viewer draws no profile layers
    for that structure.
    """
    def q(p):
        return np.round((np.asarray(p, dtype=np.float64) - center) * POS_SCALE).astype(np.int16)

    out = {
        'nframes': 1,
        'x1': {'counts': b64(np.array([len(el)], np.uint16)),
               'elem': b64(np.array([pal.elem(e) for e in el], np.uint8)),
               'pos': b64(q(xyz))},
    }
    if bonds:
        out['x1']['bonds'] = b64(np.array(bonds, np.uint16).reshape(-1, 2))

    if s3 is None:
        return out

    pos3 = np.array([s[1] for s in s3]) if s3 else np.zeros((0, 3))
    esp = np.array([s[0][2] - s[0][0] for s in s3]) if s3 else np.zeros(0)
    out['x3'] = {'n': len(s3),
                 'pos': b64(q(pos3)),
                 'esp': b64(np.round(np.clip(esp, -1, 1) * 127).astype(np.int8))}

    out['x4s'] = {'counts': b64(np.array([len(s4)], np.uint16)),
                  'col': b64(np.array([pal.color(c) for c, _, _ in s4], np.uint8)),
                  'rad': b64(np.array([int(round(r * 50)) for _, _, r in s4], np.uint8)),
                  'pos': b64(q([p for _, p, _ in s4]).reshape(-1, 3))}

    dirs, lens = [], []
    for _, p0, p1 in a4:
        v = np.array(p1) - np.array(p0)
        L = np.linalg.norm(v) or 1.0
        dirs.append(np.round(v / L * 127).astype(np.int8))
        lens.append(int(min(255, round(L * 20))))
    out['x4a'] = {'counts': b64(np.array([len(a4)], np.uint16)),
                  'col': b64(np.array([pal.color(c) for c, _, _ in a4], np.uint8)),
                  'len': b64(np.array(lens, np.uint8)),
                  'pos': b64(q([p0 for _, p0, _ in a4]).reshape(-1, 3)),
                  'dir': b64(np.array(dirs, np.int8).reshape(-1, 3))}
    return out


def pack_trajectory(frame_dir, stride, center, pal):
    """Atoms only.

    The x3 and x4 channels do denoise alongside the atoms, but a surface and a
    set of pharmacophores that are still mostly noise say nothing a viewer can
    read, and they bury the molecule they belong to.  The profiles that matter
    are the fragments' -- those are a fixed target, and they are packed with
    the static structures.  Dropping the diffusing ones here is also what takes
    each trajectory from ~180 KB to a fraction of it.
    """
    idx = sorted(int(re.search(r'(\d+)_frame_x1\.xyz$', p).group(1))
                 for p in glob.glob(os.path.join(frame_dir, '*_frame_x1.xyz')))
    if not idx:
        sys.exit(f'no frames in {frame_dir}')
    keep = sorted(set(idx[::stride]) | {idx[-1]})   # the final frame always survives

    def q(p):
        return np.round((np.asarray(p, dtype=np.float64) - center) * POS_SCALE).astype(np.int16)

    x1_counts, x1_elem, x1_pos = [], [], []

    for i in keep:
        s = f'{i:03d}'
        el, p1, _ = read_xyz(os.path.join(frame_dir, f'{s}_frame_x1.xyz'))
        x1_counts.append(len(el))
        x1_elem.extend(pal.elem(e) for e in el)
        x1_pos.append(q(p1))

    return {
        'nframes': len(keep), 'srcFrames': len(idx), 'frameIndex': keep,
        'x1': {'counts': b64(np.array(x1_counts, np.uint16)),
               'elem': b64(np.array(x1_elem, np.uint8)),
               'pos': b64(np.concatenate(x1_pos))},
    }


# --------------------------------------------------------------------------
# protein ribbon
#
# Only what a ribbon actually needs leaves Python: one CA per residue and the
# carbonyl direction that orients the flat face.  The spline, the widening at
# helices and strands, and the mesh itself are built in the browser -- 182
# residues of control points is under 2 KB, where the finished mesh would be
# hundreds.
# --------------------------------------------------------------------------
def pack_protein(pdb, center):
    res = {}                      # (chain, resseq) -> {CA, O}
    order = []
    for l in open(pdb):
        if not l.startswith('ATOM'):
            continue
        name, chain, seq = l[12:16].strip(), l[21], int(l[22:26])
        if name not in ('CA', 'O'):
            continue
        key = (chain, seq)
        if key not in res:
            res[key] = {}
            order.append(key)
        res[key][name] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])

    # secondary structure, so strands and helices can be drawn wider than coil
    ss = {}
    for l in open(pdb):
        if l.startswith('HELIX'):
            c = l[19]
            for s in range(int(l[21:25]), int(l[33:37]) + 1):
                ss[(c, s)] = 1
        elif l.startswith('SHEET'):
            c = l[21]
            for s in range(int(l[22:26]), int(l[33:37]) + 1):
                ss[(c, s)] = 2

    keys = [k for k in order if 'CA' in res[k]]
    ca = np.array([res[k]['CA'] for k in keys])
    # carbonyl direction orients the ribbon's flat face; the last residue of a
    # run has no usable O, so it inherits its neighbour's.
    d = np.zeros_like(ca)
    for i, k in enumerate(keys):
        if 'O' in res[k]:
            v = res[k]['O'] - res[k]['CA']
        else:
            v = d[i - 1] if i else np.array([0.0, 0.0, 1.0])
        n = np.linalg.norm(v) or 1.0
        d[i] = v / n
    # Sign of the carbonyl flips residue to residue along a strand; without
    # this the ribbon would corkscrew 180 degrees between every pair.
    for i in range(1, len(d)):
        if d[i] @ d[i - 1] < 0:
            d[i] = -d[i]

    # a gap in the model is a real break -- the ribbon must not bridge it
    gap = np.zeros(len(ca), np.uint8)
    for i in range(1, len(ca)):
        if np.linalg.norm(ca[i] - ca[i - 1]) > 4.5 or keys[i][0] != keys[i - 1][0]:
            gap[i] = 1

    return {
        'n': len(ca),
        'scale': PRO_SCALE,
        'ca': b64(np.round((ca - center) * PRO_SCALE).astype(np.int16)),
        'dir': b64(np.round(d * 127).astype(np.int8)),
        'ss': b64(np.array([ss.get(k, 0) for k in keys], np.uint8)),
        'brk': b64(gap),
    }


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('merge_dir')
    ap.add_argument('--pdb', default='assets/data/8cnx_apo.pdb')
    ap.add_argument('--out', default='assets/data')
    ap.add_argument('--stride', type=int, default=3)
    a = ap.parse_args()
    D = a.merge_dir

    statics = {
        'frag8':   ('chimera_frag8', '8', 'Fragment 8'),
        'frag10':  ('chimera_frag10', '10', 'Fragment 10'),
        'frag12':  ('chimera_frag12', '12', 'Fragment 12'),
        'merged1': ('chimera_merged_8_10', 'merged_8_10', 'Merge of 8 + 10'),
    }

    # ---- one centre for the whole story: the middle of everything that binds
    pts = []
    for sub, stem, _ in statics.values():
        _, xyz, _ = read_sdf(os.path.join(D, sub, f'{stem}_x1.sdf'))
        pts.append(xyz)
    r2_el, r2_xyz, _ = read_xyz(os.path.join(D, 'sample_round_2', '359_frame_x1.xyz'))
    pts.append(r2_xyz)
    center = np.concatenate(pts).mean(axis=0)
    print('pocket centre:', center.round(3))

    pal = Palette()
    out = {'center': [round(float(v), 4) for v in center],
           'posScale': POS_SCALE,
           'protein': pack_protein(a.pdb, center),
           'structures': {}}

    for key, (sub, stem, label) in statics.items():
        el, xyz, bonds = read_sdf(os.path.join(D, sub, f'{stem}_x1.sdf'))
        s3, _ = read_bild(os.path.join(D, sub, f'{stem}_x3.bild'))
        s4, a4 = read_bild(os.path.join(D, sub, f'{stem}_x4.bild'))
        out['structures'][key] = pack_structure(el, xyz, bonds, s3, s4, a4, center, pal)
        print(f'  {key:8s} {len(el):3d} atoms, {len(s3)} ESP pts, '
              f'{len(s4)} pharmacophores, {len(a4)} vectors')

    # Round 2 has no .sdf -- its product is the last frame of its own
    # trajectory, and so is the only profile it could carry: the last step of a
    # denoising run rather than a profile computed from a finished molecule.
    # It is not the same object as a fragment's measured profile, so it is left
    # out and the final stage shows the molecule alone.
    out['structures']['merged2'] = pack_structure(r2_el, r2_xyz, None, None, None, None,
                                                  center, pal)
    print(f'  merged2  {len(r2_el):3d} atoms, no profile')

    trajs = {
        'evo-traj-1.js': ('EVO_TRAJ_1', os.path.join(D, 'sample_0', 'x')),
        'evo-traj-2.js': ('EVO_TRAJ_2', os.path.join(D, 'sample_round_2')),
    }
    packed = {}
    for fn, (var, path) in trajs.items():
        packed[fn] = (var, pack_trajectory(path, a.stride, center, pal))
        print(f'  {fn}: {packed[fn][1]["nframes"]} of '
              f'{packed[fn][1]["srcFrames"]} frames (stride {a.stride})')

    # the palettes are shared, so they can only be written once every structure
    # and every trajectory frame has been through them
    out['elements'] = pal.elements
    out['colors'] = [list(c) for c in pal.colors]
    for fn, (var, p) in packed.items():
        p['posScale'] = POS_SCALE
        p['elements'] = pal.elements
        p['colors'] = [list(c) for c in pal.colors]

    def write(fn, var, payload):
        path = os.path.join(a.out, fn)
        with open(path, 'w') as fh:
            fh.write('/* generated by tools/build_evolution.py -- do not edit */\n')
            fh.write(f'window.{var} = ')
            json.dump(payload, fh, separators=(',', ':'))
            fh.write(';\n')
        print(f'wrote {path}  ({os.path.getsize(path)/1024:.0f} KB)')

    write('evo-base.js', 'EVO_BASE', out)
    for fn, (var, p) in packed.items():
        write(fn, var, p)
    print('elements:', pal.elements)
    print('colours :', pal.colors)


if __name__ == '__main__':
    main()
