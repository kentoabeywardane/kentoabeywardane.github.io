/* Shared molecule rendering ---------------------------------------------------
   One "actor" is one thing the scene can show: a crystallographic fragment, a
   generated molecule, or a whole denoising trajectory.  Statics and
   trajectories use the same payload layout -- a static is simply a trajectory
   with one frame -- so a single draw path covers both and the staged viewer
   can treat every object in the scene alike.

   Each actor owns three layers:
     x1  the molecule   (atom spheres + inferred or supplied bonds)
     x3  the ESP surface(shell of points coloured by signed potential)
     x4  the pharmacophores (typed anchors + direction vectors)
---------------------------------------------------------------------------- */
import * as THREE from 'three';

/* CPK, slightly desaturated so several molecules can sit on top of each other
   without the scene turning into a colour fight. */
export const ELEM = {
  H:  { c: 0xe8e8ea, r: 0.30 }, C:  { c: 0x4c5057, r: 0.46 },
  N:  { c: 0x3b6fd4, r: 0.44 }, O:  { c: 0xd4453b, r: 0.43 },
  F:  { c: 0x5bbf7a, r: 0.40 }, P:  { c: 0xe08a2e, r: 0.52 },
  S:  { c: 0xd9b63a, r: 0.52 }, Cl: { c: 0x4fae5a, r: 0.50 },
  Br: { c: 0xa2543a, r: 0.55 }, I:  { c: 0x8c4fc4, r: 0.60 },
  Si: { c: 0xb0a090, r: 0.56 },
};
export const FALLBACK = { c: 0x9aa0a8, r: 0.45 };
const COV = { H: 0.31, C: 0.76, N: 0.71, O: 0.66, F: 0.57, P: 1.07, S: 1.05,
              Cl: 1.02, Br: 1.20, I: 1.39, Si: 1.11 };

export function unb64(s, Type) {
  const bin = atob(s);
  const buf = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) buf[i] = bin.charCodeAt(i);
  return new Type(buf.buffer, 0, bin.length / Type.BYTES_PER_ELEMENT);
}
function offsets(counts) {
  const o = new Uint32Array(counts.length + 1);
  for (let i = 0; i < counts.length; i++) o[i + 1] = o[i] + counts[i];
  return o;
}
/* diverging ESP ramp: red (negative) -> neutral -> blue (positive) */
export function espColor(t, out) {
  const a = Math.abs(t);
  const lo = [0.88, 0.88, 0.86];
  const hi = t >= 0 ? [0.11, 0.38, 0.85] : [0.84, 0.17, 0.14];
  out.setRGB(lo[0] + (hi[0] - lo[0]) * a,
             lo[1] + (hi[1] - lo[1]) * a,
             lo[2] + (hi[2] - lo[2]) * a, THREE.SRGBColorSpace);
  return out;
}

/* ---------------------------------------------------------------------------
   createActor(data, opts)

   opts.carbon    hex that replaces carbon's colour, so two molecules sharing a
                  scene can be told apart the way a medicinal-chemistry figure
                  does it -- heteroatoms keep their CPK colours either way.
   opts.scale     radius multiplier; overlaid references are drawn thinner than
                  the molecule the stage is actually about.
   opts.palette / opts.elements  shared tables when the payload does not carry
                  its own (the evolution payloads share one set across files).
--------------------------------------------------------------------------- */
export function createActor(data, opts = {}) {
  const S = data.posScale;
  const nF = data.nframes;
  const scale = opts.scale ?? 1;

  const x1c = unb64(data.x1.counts, Uint16Array), x1o = offsets(x1c);
  const x1e = unb64(data.x1.elem, Uint8Array);
  const x1p = unb64(data.x1.pos, Int16Array);
  const fixedBonds = data.x1.bonds ? unb64(data.x1.bonds, Uint16Array) : null;

  const n3 = data.x3.n;
  const x3p = unb64(data.x3.pos, Int16Array);
  const x3e = unb64(data.x3.esp, Int8Array);

  const s4c = unb64(data.x4s.counts, Uint16Array), s4o = offsets(s4c);
  const s4col = unb64(data.x4s.col, Uint8Array);
  const s4rad = unb64(data.x4s.rad, Uint8Array);
  const s4p = unb64(data.x4s.pos, Int16Array);

  const a4c = unb64(data.x4a.counts, Uint16Array), a4o = offsets(a4c);
  const a4col = unb64(data.x4a.col, Uint8Array);
  const a4len = unb64(data.x4a.len, Uint8Array);
  const a4p = unb64(data.x4a.pos, Int16Array);
  const a4d = unb64(data.x4a.dir, Int8Array);

  const elementNames = opts.elements || data.elements;
  const palette = (opts.colors || data.colors).map(
    c => new THREE.Color().setRGB(c[0], c[1], c[2], THREE.SRGBColorSpace));
  const elemStyle = elementNames.map((e) => {
    const base = ELEM[e] || FALLBACK;
    return (e === 'C' && opts.carbon != null) ? { c: opts.carbon, r: base.r } : base;
  });
  const elemCov = elementNames.map(e => COV[e] ?? 0.9);
  const isH = elementNames.map(e => e === 'H');

  const MAXA = Math.max(...x1c);
  const MAXS = Math.max(1, ...s4c);
  const MAXR = Math.max(1, ...a4c);
  const MAXB = fixedBonds ? fixedBonds.length : MAXA * 4;

  const mkMat = (o) => new THREE.MeshStandardMaterial(
    Object.assign({ roughness: 0.42, metalness: 0.03 }, o));

  const atoms = new THREE.InstancedMesh(
    new THREE.SphereGeometry(1, 20, 14), mkMat({}), MAXA);
  const bonds = new THREE.InstancedMesh(
    new THREE.CylinderGeometry(0.11, 0.11, 1, 10), mkMat({}), MAXB);
  const esp = new THREE.InstancedMesh(
    new THREE.SphereGeometry(1, 12, 9),
    mkMat({ transparent: true, opacity: 0.52, roughness: 0.6 }), Math.max(1, n3));
  const pharm = new THREE.InstancedMesh(
    new THREE.SphereGeometry(1, 18, 12),
    mkMat({ transparent: true, opacity: 0.4, depthWrite: false }), MAXS);
  const shafts = new THREE.InstancedMesh(
    new THREE.CylinderGeometry(0.1, 0.1, 1, 8), mkMat({}), MAXR);
  const heads = new THREE.InstancedMesh(
    new THREE.ConeGeometry(0.2, 0.46, 10), mkMat({}), MAXR);

  const meshes = [atoms, bonds, esp, pharm, shafts, heads];
  for (const m of meshes) m.frustumCulled = false;

  const group = new THREE.Group();
  const layers = {
    x1: new THREE.Group(), x3: new THREE.Group(), x4: new THREE.Group(),
  };
  layers.x1.add(atoms, bonds);
  layers.x3.add(esp);
  layers.x4.add(pharm, shafts, heads);
  group.add(layers.x1, layers.x3, layers.x4);

  /* ---------- drawing ---------- */
  const M = new THREE.Matrix4(), Q = new THREE.Quaternion();
  const A = new THREE.Vector3(), B = new THREE.Vector3(), V = new THREE.Vector3();
  const UP = new THREE.Vector3(0, 1, 0);
  const C = new THREE.Color();
  const px = new Float32Array(MAXA * 3);
  const keep = new Int32Array(MAXA);        /* payload index -> draw slot, or -1 */

  function placeCylinder(mesh, i, x0, y0, z0, x1_, y1_, z1_, radius) {
    A.set(x0, y0, z0); B.set(x1_, y1_, z1_);
    V.subVectors(B, A);
    const len = V.length() || 1e-6;
    Q.setFromUnitVectors(UP, V.divideScalar(len));
    M.compose(A.addVectors(A, V.multiplyScalar(len / 2)), Q,
              new THREE.Vector3(radius, len, radius));
    mesh.setMatrixAt(i, M);
  }

  let frame = -1, showH = opts.showH ?? false;

  function draw(f) {
    frame = f;
    /* --- atoms --- */
    const a0 = x1o[f], na = x1c[f];
    let k = 0;
    for (let i = 0; i < na; i++) {
      const e = x1e[a0 + i];
      if (!showH && isH[e]) { keep[i] = -1; continue; }
      const st = elemStyle[e];
      const x = x1p[(a0 + i) * 3] / S, y = x1p[(a0 + i) * 3 + 1] / S,
            z = x1p[(a0 + i) * 3 + 2] / S;
      px[k * 3] = x; px[k * 3 + 1] = y; px[k * 3 + 2] = z;
      M.compose(A.set(x, y, z), Q.identity(), V.setScalar(st.r * 0.55 * scale));
      atoms.setMatrixAt(k, M);
      atoms.setColorAt(k, C.setHex(st.c));
      keep[i] = k++;
    }
    atoms.count = k;
    atoms.instanceMatrix.needsUpdate = true;
    if (atoms.instanceColor) atoms.instanceColor.needsUpdate = true;

    /* --- bonds --- */
    let nb = 0;
    const halfBond = (xi, yi, zi, xj, yj, zj, ci, cj) => {
      const mx = (xi + xj) / 2, my = (yi + yj) / 2, mz = (zi + zj) / 2;
      placeCylinder(bonds, nb, xi, yi, zi, mx, my, mz, scale);
      bonds.setColorAt(nb++, C.setHex(ci));
      placeCylinder(bonds, nb, mx, my, mz, xj, yj, zj, scale);
      bonds.setColorAt(nb++, C.setHex(cj));
    };

    if (fixedBonds) {
      /* a real bond table from the .sdf: no guessing needed */
      for (let b = 0; b < fixedBonds.length / 2 && nb < MAXB - 1; b++) {
        const i = fixedBonds[b * 2], j = fixedBonds[b * 2 + 1];
        const ki = keep[i], kj = keep[j];
        if (ki < 0 || kj < 0) continue;
        halfBond(px[ki * 3], px[ki * 3 + 1], px[ki * 3 + 2],
                 px[kj * 3], px[kj * 3 + 1], px[kj * 3 + 2],
                 elemStyle[x1e[a0 + i]].c, elemStyle[x1e[a0 + j]].c);
      }
    } else {
      /* trajectory frames carry no connectivity, so bonds come from a
         distance test against covalent radii -- which is why they flicker
         early on, while the geometry is still noise */
      const el = [];
      for (let i = 0; i < na; i++) if (keep[i] >= 0) el.push(x1e[a0 + i]);
      for (let i = 0; i < k && nb < MAXB - 1; i++) {
        const ri = elemCov[el[i]];
        for (let j = i + 1; j < k && nb < MAXB - 1; j++) {
          const dx = px[i * 3] - px[j * 3], dy = px[i * 3 + 1] - px[j * 3 + 1],
                dz = px[i * 3 + 2] - px[j * 3 + 2];
          const d2 = dx * dx + dy * dy + dz * dz;
          const cut = ri + elemCov[el[j]] + 0.4;
          if (d2 > cut * cut || d2 < 0.16) continue;
          halfBond(px[i * 3], px[i * 3 + 1], px[i * 3 + 2],
                   px[j * 3], px[j * 3 + 1], px[j * 3 + 2],
                   elemStyle[el[i]].c, elemStyle[el[j]].c);
        }
      }
    }
    bonds.count = nb;
    bonds.instanceMatrix.needsUpdate = true;
    if (bonds.instanceColor) bonds.instanceColor.needsUpdate = true;

    /* --- ESP surface --- */
    const e0 = f * n3;
    for (let i = 0; i < n3; i++) {
      M.compose(A.set(x3p[(e0 + i) * 3] / S, x3p[(e0 + i) * 3 + 1] / S,
                      x3p[(e0 + i) * 3 + 2] / S),
                Q.identity(), V.setScalar(0.15));
      esp.setMatrixAt(i, M);
      esp.setColorAt(i, espColor(x3e[e0 + i] / 127, C));
    }
    esp.count = n3;
    esp.instanceMatrix.needsUpdate = true;
    if (esp.instanceColor) esp.instanceColor.needsUpdate = true;

    /* --- pharmacophore anchors --- */
    const p0 = s4o[f], ns = s4c[f];
    for (let i = 0; i < ns; i++) {
      M.compose(A.set(s4p[(p0 + i) * 3] / S, s4p[(p0 + i) * 3 + 1] / S,
                      s4p[(p0 + i) * 3 + 2] / S),
                Q.identity(), V.setScalar(Math.max(0.3, s4rad[p0 + i] / 50)));
      pharm.setMatrixAt(i, M);
      pharm.setColorAt(i, palette[s4col[p0 + i]]);
    }
    pharm.count = ns;
    pharm.instanceMatrix.needsUpdate = true;
    if (pharm.instanceColor) pharm.instanceColor.needsUpdate = true;

    /* --- pharmacophore direction vectors --- */
    const r0 = a4o[f], nr = a4c[f];
    for (let i = 0; i < nr; i++) {
      const j = r0 + i;
      const ox = a4p[j * 3] / S, oy = a4p[j * 3 + 1] / S, oz = a4p[j * 3 + 2] / S;
      const L = a4len[j] / 20;
      const dx = a4d[j * 3] / 127, dy = a4d[j * 3 + 1] / 127, dz = a4d[j * 3 + 2] / 127;
      const sl = Math.max(0.01, L - 0.46);
      placeCylinder(shafts, i, ox, oy, oz, ox + dx * sl, oy + dy * sl, oz + dz * sl, 1);
      shafts.setColorAt(i, palette[a4col[j]]);
      V.set(dx, dy, dz).normalize();
      Q.setFromUnitVectors(UP, V);
      M.compose(A.set(ox + dx * (L - 0.23), oy + dy * (L - 0.23), oz + dz * (L - 0.23)),
                Q, B.set(1, 1, 1));
      heads.setMatrixAt(i, M);
      heads.setColorAt(i, palette[a4col[j]]);
    }
    shafts.count = nr; heads.count = nr;
    shafts.instanceMatrix.needsUpdate = true; heads.instanceMatrix.needsUpdate = true;
    if (shafts.instanceColor) shafts.instanceColor.needsUpdate = true;
    if (heads.instanceColor) heads.instanceColor.needsUpdate = true;
  }

  /* ---------- fades ----------
     Opacity is driven per layer so a stage can bring a reference molecule up
     as a ghost without it popping into existence. */
  const baseOpacity = new Map([
    [atoms, 1], [bonds, 1], [esp, 0.52], [pharm, 0.4], [shafts, 1], [heads, 1],
  ]);
  function setLayerOpacity(layer, k) {
    for (const m of layers[layer].children) {
      const base = baseOpacity.get(m) * k;
      m.material.opacity = base;
      m.material.transparent = base < 0.999;
      m.material.depthWrite = m === pharm ? false : base > 0.6;
      m.visible = base > 0.004;
    }
  }

  /* bounding radius about the world origin, for camera framing */
  function radius(f = Math.max(0, frame), withSurface = true) {
    let r = 0;
    const a0 = x1o[f];
    for (let i = 0; i < x1c[f]; i++)
      r = Math.max(r, Math.hypot(x1p[(a0 + i) * 3] / S, x1p[(a0 + i) * 3 + 1] / S,
                                 x1p[(a0 + i) * 3 + 2] / S));
    if (withSurface) {
      const e0 = f * n3;
      for (let i = 0; i < n3; i++)
        r = Math.max(r, Math.hypot(x3p[(e0 + i) * 3] / S, x3p[(e0 + i) * 3 + 1] / S,
                                   x3p[(e0 + i) * 3 + 2] / S));
    }
    return r;
  }

  function dispose() {
    for (const m of meshes) { m.geometry.dispose(); m.material.dispose(); m.dispose(); }
  }

  draw(0);
  return {
    group, layers, nframes: nF, data,
    setFrame: (f) => draw(Math.max(0, Math.min(nF - 1, f | 0))),
    setShowH: (v) => { if (v !== showH) { showH = v; draw(Math.max(0, frame)); } },
    setLayerOpacity, radius, dispose,
    get frame() { return frame; },
  };
}
