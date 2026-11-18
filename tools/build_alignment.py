#!/usr/bin/env python3
"""Pack a shepherd_score alignment-trajectory dump into a JS payload for the
alignment viewer (shepherd-score-alignment.html).

Usage (needs rdkit + shepherd_score, e.g. the shepherdscore conda env):
  python tools/build_alignment.py <pair_dir> <molblock_pickle> <out.js> [--var NAME]

<pair_dir> is a ``mol<i>__mol<j>/`` folder written by gen_alignment_trajs.py:
  info.json, ref_x1.xyz, ref_x3.bild
  vol_noH/NNN_frame_x1.xyz (heavy atoms)    scores.json
  vol_avoid/NNN_frame_x1.xyz (heavy atoms)  scores.json (+ similarity, penalty)
  avoid_points.xyz                          (vol_avoid exclusion wall; add_vol_avoid.py)
  surf/NNN_frame_x1.xyz, NNN_frame_x3.bild  scores.json
  pharm/NNN_frame_x1.xyz                    scores.json

The fit molecule moves rigidly, so instead of shipping every frame's coordinates
the payload keeps the fit once (in its start pose) and, per frame, the rigid
transform from the start pose (quaternion xyzw + translation), recovered by
Kabsch from the dumped coordinates. Pharmacophores are recomputed with
shepherd_score from the same molblocks (deterministic); surface points are read
from the BILD files because sampling them is stochastic.
"""
import argparse, glob, json, os, pickle, re
import numpy as np
from rdkit import Chem
from scipy.spatial.transform import Rotation
from shepherd_score.container import Molecule

MODES = ('vol_noH', 'vol_avoid', 'surf', 'pharm')
HEAVY_ONLY = ('vol_noH', 'vol_avoid')


def read_xyz(path):
    lines = [l for l in open(path).read().splitlines() if l.strip()]
    n = int(lines[0])
    rows = lines[1:][-n:]
    return [r.split()[0] for r in rows], np.array([[float(v) for v in r.split()[1:4]] for r in rows])


def read_spheres(path):
    return np.array([[float(v) for v in l.split()[1:4]] for l in open(path) if l.startswith('.sphere')])


def kabsch(A, B):
    """R, t minimizing |R a + t - b|."""
    ca, cb = A.mean(0), B.mean(0)
    U, _, Vt = np.linalg.svd((A - ca).T @ (B - cb))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, cb - R @ ca


def frames(d):
    return sorted(glob.glob(os.path.join(d, '*_frame_x1.xyz')),
                  key=lambda p: int(re.search(r'(\d+)_frame', os.path.basename(p)).group(1)))


def r(a, n=3):
    return np.round(np.asarray(a, dtype=float), n).tolist()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('pair_dir')
    ap.add_argument('pickle')
    ap.add_argument('out')
    ap.add_argument('--var', default='SHEPHERD_ALIGN')
    a = ap.parse_args()

    info = json.load(open(os.path.join(a.pair_dir, 'info.json')))
    data = pickle.load(open(a.pickle, 'rb'))

    def mol_and_pharm(idx, coords):
        mol = Chem.MolFromMolBlock(data[idx][0], removeHs=False)
        conf = mol.GetConformer()
        for k, p in enumerate(coords):
            conf.SetAtomPosition(k, p.tolist())
        m = Molecule(mol, pharm_multi_vector=False)
        return mol, m

    # reference, as dumped (already centered)
    _, ref_pos = read_xyz(os.path.join(a.pair_dir, 'ref_x1.xyz'))
    ref_mol, ref_m = mol_and_pharm(info['ref_index'], ref_pos)
    # fit, in the shared start pose (frame 0 of the all-atom tracks)
    _, fit_pos = read_xyz(frames(os.path.join(a.pair_dir, 'surf'))[0])
    fit_mol, fit_m = mol_and_pharm(info['fit_index'], fit_pos)
    heavy = np.array([at.GetAtomicNum() > 1 for at in fit_mol.GetAtoms()])

    def mol_json(mol, pos, m, surf):
        return dict(
            elements=[at.GetSymbol() for at in mol.GetAtoms()],
            pos=r(pos),
            bonds=[[b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondTypeAsDouble()] for b in mol.GetBonds()],
            surf=r(surf),
            pharm=dict(types=[int(t) for t in m.pharm_types], ancs=r(m.pharm_ancs), vecs=r(m.pharm_vecs)),
        )

    ref_surf = read_spheres(os.path.join(a.pair_dir, 'ref_x3.bild'))
    fit_surf = read_spheres(os.path.join(a.pair_dir, 'surf', '000_frame_x3.bild'))

    modes = {}
    for mode in MODES:
        d = os.path.join(a.pair_dir, mode)
        start = fit_pos[heavy] if mode in HEAVY_ONLY else fit_pos
        quat, trans, err = [], [], 0.0
        for p in frames(d):
            _, X = read_xyz(p)
            R, t = kabsch(start, X)
            err = max(err, float(np.abs(start @ R.T + t - X).max()))
            quat.append(r(Rotation.from_matrix(R).as_quat(), 5))
            trans.append(r(t, 4))
        sj = json.load(open(os.path.join(d, 'scores.json')))
        scores = sj['scores']
        assert len(scores) == len(quat), (mode, len(scores), len(quat))
        modes[mode] = dict(quat=quat, trans=trans, scores=r(scores, 4))
        if 'penalty' in sj:
            modes[mode].update(similarity=r(sj['similarity'], 4), penalty=r(sj['penalty'], 4))
        print(f'{mode:8s} {len(quat):3d} frames  score {scores[0]:.3f} -> {scores[-1]:.3f}  max Kabsch residual {err:.1e} A')

    out = dict(
        info={k: info[k] for k in ('ref_index', 'fit_index', 'ref_smiles', 'fit_smiles',
                                   'start_rotation_deg', 'start_shift_A', 'library_best')},
        ref=mol_json(ref_mol, ref_pos, ref_m, ref_surf),
        fit=mol_json(fit_mol, fit_pos, fit_m, fit_surf),
        modes=modes,
        avoid=r(read_xyz(os.path.join(a.pair_dir, 'avoid_points.xyz'))[1]),
    )
    js = f'window.{a.var} = {json.dumps(out, separators=(",", ":"))};\n'
    open(a.out, 'w').write(js)
    print(f'wrote {a.out} ({len(js) / 1024:.0f} KB)')


if __name__ == '__main__':
    main()
