/* Protein cartoon ribbon ------------------------------------------------------
   Builds a single BufferGeometry from the compact CA trace written by
   tools/build_evolution.py.  The mesh is a rectangular cross-section swept
   along a Catmull-Rom spline through the alpha carbons: wide and flat through
   helices and strands, narrow and nearly round through coil, so the fold reads
   at a glance without anything as heavy as a real cartoon renderer.

   The flat face is oriented by the carbonyl direction stored per residue,
   which is what stops a ribbon from corkscrewing along a strand.
---------------------------------------------------------------------------- */
import * as THREE from 'three';

const SUB = 8;          /* spline samples per residue */
const HW = [0.30, 1.05, 1.15];   /* half-width: coil, helix, strand */
const HT = [0.30, 0.22, 0.22];   /* half-thickness, same order */

function unb64(s, Type) {
  const bin = atob(s);
  const buf = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) buf[i] = bin.charCodeAt(i);
  return new Type(buf.buffer, 0, bin.length / Type.BYTES_PER_ELEMENT);
}

/* Everything the viewer shows sits in one pocket, and the rest of the protein
   is context.  Rather than clipping it off, the ribbon fades with distance
   from the world origin: solid around the binding site, nearly gone at the far
   side of the fold, so the molecules are never read through five stacked
   layers of translucent backbone. */
export function buildRibbon(p, o = {}) {
  const S = p.scale;
  /* only the alpha channel varies per vertex; the tint is the material's
     colour, which keeps this out of the colour-space question entirely */
  const nearR = o.nearR ?? 11, farR = o.farR ?? 30;
  const nearA = o.nearA ?? 1.0, farA = o.farA ?? 0.34;
  const alphaAt = (x, y, z) => {
    const t = THREE.MathUtils.clamp((Math.hypot(x, y, z) - nearR) / (farR - nearR), 0, 1);
    return nearA + (farA - nearA) * (t * t * (3 - 2 * t));   /* smoothstep */
  };
  const ca = unb64(p.ca, Int16Array);
  const dv = unb64(p.dir, Int8Array);
  const ss = unb64(p.ss, Uint8Array);
  const brk = unb64(p.brk, Uint8Array);
  const n = p.n;

  /* split into unbroken runs -- a gap in the model is not a bond */
  const runs = [];
  let start = 0;
  for (let i = 1; i <= n; i++) {
    if (i === n || brk[i]) {
      if (i - start >= 2) runs.push([start, i]);
      start = i;
    }
  }

  const pos = [], nrm = [], col = [], idx = [];
  const P = new THREE.Vector3(), T = new THREE.Vector3(), W = new THREE.Vector3(),
        N = new THREE.Vector3(), tmp = new THREE.Vector3();

  for (const [a, b] of runs) {
    const ptsA = [], ptsB = [];
    for (let i = a; i < b; i++) {
      const x = ca[i * 3] / S, y = ca[i * 3 + 1] / S, z = ca[i * 3 + 2] / S;
      ptsA.push(new THREE.Vector3(x, y, z));
      /* a second curve offset along the carbonyl: interpolating both and
         subtracting keeps the orientation as smooth as the path itself */
      ptsB.push(new THREE.Vector3(x + dv[i * 3] / 127, y + dv[i * 3 + 1] / 127,
                                  z + dv[i * 3 + 2] / 127));
    }
    const cA = new THREE.CatmullRomCurve3(ptsA, false, 'catmullrom', 0.5);
    const cB = new THREE.CatmullRomCurve3(ptsB, false, 'catmullrom', 0.5);

    const nres = b - a;
    const steps = (nres - 1) * SUB;
    const ring0 = pos.length / 3;

    for (let s = 0; s <= steps; s++) {
      const t = s / steps;
      cA.getPoint(t, P);
      cA.getTangent(t, T).normalize();
      W.copy(cB.getPoint(t, tmp)).sub(P);
      /* make the face direction square to the path, then complete the frame */
      W.addScaledVector(T, -W.dot(T));
      if (W.lengthSq() < 1e-8) W.set(T.y, -T.x, 0);      /* degenerate guard */
      W.normalize();
      N.crossVectors(T, W).normalize();

      /* width follows secondary structure, smoothed between residues */
      const g = t * (nres - 1), i0 = Math.min(nres - 1, Math.floor(g));
      const i1 = Math.min(nres - 1, i0 + 1), fr = g - i0;
      const k0 = ss[a + i0], k1 = ss[a + i1];
      const hw = HW[k0] + (HW[k1] - HW[k0]) * fr;
      const ht = HT[k0] + (HT[k1] - HT[k0]) * fr;
      /* taper the open ends so a run finishes in a point rather than a stub */
      const taper = Math.min(1, Math.min(s, steps - s) / (SUB * 0.6) * 0.75 + 0.25);

      const alpha = alphaAt(P.x, P.y, P.z);
      for (const [sw, sn] of [[1, 1], [1, -1], [-1, -1], [-1, 1]]) {
        pos.push(P.x + W.x * hw * sw * taper + N.x * ht * sn,
                 P.y + W.y * hw * sw * taper + N.y * ht * sn,
                 P.z + W.z * hw * sw * taper + N.z * ht * sn);
        col.push(1, 1, 1, alpha);
      }
      /* corner normals: the average of the two faces that meet there, which
         rounds the edge just enough to catch the light */
      for (const [sw, sn] of [[1, 1], [1, -1], [-1, -1], [-1, 1]]) {
        tmp.copy(W).multiplyScalar(sw).addScaledVector(N, sn).normalize();
        nrm.push(tmp.x, tmp.y, tmp.z);
      }
    }

    for (let s = 0; s < steps; s++) {
      const r0 = ring0 + s * 4, r1 = r0 + 4;
      for (let c = 0; c < 4; c++) {
        const c1 = (c + 1) % 4;
        idx.push(r0 + c, r1 + c, r0 + c1, r0 + c1, r1 + c, r1 + c1);
      }
    }
    /* caps, so an open end is not a hole you can see into */
    const last = ring0 + steps * 4;
    idx.push(ring0, ring0 + 2, ring0 + 1, ring0, ring0 + 3, ring0 + 2);
    idx.push(last, last + 1, last + 2, last, last + 2, last + 3);
  }

  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
  g.setAttribute('normal', new THREE.Float32BufferAttribute(nrm, 3));
  /* 4 components, so the alpha reaches the material -- needs vertexColors and
     transparent on a plain (non-instanced) mesh */
  g.setAttribute('color', new THREE.Float32BufferAttribute(col, 4));
  g.setIndex(idx);
  g.computeBoundingSphere();
  return g;
}
