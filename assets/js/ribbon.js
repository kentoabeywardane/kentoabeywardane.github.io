/* Protein cartoon ribbon ------------------------------------------------------
   The classical cartoon construction, as PyMOL, Chimera and NGL build it:

     1. spline the alpha-carbon trace                      (Catmull-Rom)
     2. carry a frame along it whose flat face follows the carbonyl
     3. sweep a closed cross-section through that frame
     4. widen the cross-section through helices and strands

   Three details separate this from a naive sweep, and all three are what a
   hand-rolled version gets wrong:

   * The guide direction is smoothed before it is splined.  The carbonyl of a
     single coil residue can point almost anywhere; splining the raw vectors
     makes the ribbon twist a few degrees one way and back again every
     residue, which is exactly the crimping you see along loops.  Three
     Laplacian passes along the chain take that out without touching the
     steady alternation through a strand.

   * The cross-section is a superellipse rather than a rectangle, so the
     silhouette has no corners to alias against and its normals are analytic
     and continuous.  A four-corner box has to fake its normals, and the fake
     is visible as a hard line down the middle of every face.

   * The tube is closed and wound outward, so it renders FrontSide.  A
     translucent DoubleSide tube blends its own back wall over its front wall
     and every crossing of the chain reads twice as dark as it should.

   Transparency is per vertex in the alpha channel only: solid around the
   binding site, fading with distance so the far side of the fold stays
   context instead of becoming a wall of stacked glass.  The tint is the
   material's colour, which keeps the geometry out of the colour-space
   question entirely.
---------------------------------------------------------------------------- */
import * as THREE from 'three';

const SUB = 12;                   /* spline samples per residue */
const PROFILE = 16;               /* points around the cross-section */
const SE = 5;                     /* superellipse exponent: 2 is an ellipse,
                                     large is a rectangle, 5 is a cartoon */
const HW = [0.34, 1.10, 1.25];    /* half-width:     coil, helix, strand */
const HT = [0.34, 0.26, 0.26];    /* half-thickness, same order */
const SMOOTH_PASSES = 3;

function unb64(s, Type) {
  const bin = atob(s);
  const buf = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) buf[i] = bin.charCodeAt(i);
  return new Type(buf.buffer, 0, bin.length / Type.BYTES_PER_ELEMENT);
}

/* One closed cross-section, as offsets and outward normals in the (width,
   thickness) plane of the frame.  Superellipse |u/a|^n + |v/b|^n = 1, whose
   gradient gives the normal in closed form. */
function sectionProfile() {
  const uu = [], vv = [], nu = [], nv = [];
  for (let k = 0; k < PROFILE; k++) {
    const t = (k / PROFILE) * Math.PI * 2;
    const c = Math.cos(t), s = Math.sin(t);
    const u = Math.sign(c) * Math.pow(Math.abs(c), 2 / SE);
    const v = Math.sign(s) * Math.pow(Math.abs(s), 2 / SE);
    uu.push(u); vv.push(v);
    /* gradient of the implicit form, before the per-sample a/b scaling */
    nu.push(Math.sign(u) * Math.pow(Math.abs(u), SE - 1));
    nv.push(Math.sign(v) * Math.pow(Math.abs(v), SE - 1));
  }
  return { uu, vv, nu, nv };
}
const SEC = sectionProfile();

const smoothstep = (t) => t * t * (3 - 2 * t);

export function buildRibbon(p, o = {}) {
  const S = p.scale;
  const nearR = o.nearR ?? 11, farR = o.farR ?? 30;
  const nearA = o.nearA ?? 1.0, farA = o.farA ?? 0.34;
  const alphaAt = (x, y, z) => {
    const t = THREE.MathUtils.clamp((Math.hypot(x, y, z) - nearR) / (farR - nearR), 0, 1);
    return nearA + (farA - nearA) * smoothstep(t);
  };

  const ca = unb64(p.ca, Int16Array);
  const dv = unb64(p.dir, Int8Array);
  const ss = unb64(p.ss, Uint8Array);
  const brk = unb64(p.brk, Uint8Array);
  const n = p.n;

  /* unbroken runs: a gap in the model is not a bond */
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
        N = new THREE.Vector3(), tmp = new THREE.Vector3(), prevW = new THREE.Vector3();
  const triA = new THREE.Vector3(), triB = new THREE.Vector3(), triN = new THREE.Vector3();

  for (const [a, b] of runs) {
    const nres = b - a;

    /* --- guide directions, smoothed along the chain ---------------------- */
    const dirs = [];
    for (let i = a; i < b; i++)
      dirs.push(new THREE.Vector3(dv[i * 3] / 127, dv[i * 3 + 1] / 127, dv[i * 3 + 2] / 127)
        .normalize());
    for (let pass = 0; pass < SMOOTH_PASSES; pass++) {
      const next = dirs.map((d, i) => {
        const acc = d.clone().multiplyScalar(2);
        if (i > 0) acc.add(dirs[i - 1]);
        if (i < dirs.length - 1) acc.add(dirs[i + 1]);
        return acc.lengthSq() < 1e-10 ? d.clone() : acc.normalize();
      });
      for (let i = 0; i < dirs.length; i++) dirs[i].copy(next[i]);
    }

    /* --- two splines: the path, and the path offset along the guide ------ */
    const ptsA = [], ptsB = [];
    for (let i = a; i < b; i++) {
      const v = new THREE.Vector3(ca[i * 3] / S, ca[i * 3 + 1] / S, ca[i * 3 + 2] / S);
      ptsA.push(v);
      ptsB.push(v.clone().add(dirs[i - a]));
    }
    const cA = new THREE.CatmullRomCurve3(ptsA, false, 'catmullrom', 0.5);
    const cB = new THREE.CatmullRomCurve3(ptsB, false, 'catmullrom', 0.5);

    const steps = (nres - 1) * SUB;
    const ring0 = pos.length / 3;
    let haveWinding = false, flip = false;
    prevW.set(0, 0, 0);

    for (let s = 0; s <= steps; s++) {
      const t = s / steps;
      cA.getPoint(t, P);
      cA.getTangent(t, T).normalize();
      W.copy(cB.getPoint(t, tmp)).sub(P);
      W.addScaledVector(T, -W.dot(T));                 /* square to the path */
      if (W.lengthSq() < 1e-8) W.set(T.y, -T.x, 0);    /* degenerate guard */
      W.normalize();
      /* never let the face swap sides between two adjacent samples */
      if (prevW.lengthSq() > 0 && W.dot(prevW) < 0) W.negate();
      prevW.copy(W);
      N.crossVectors(T, W).normalize();

      /* width follows secondary structure, eased across the boundary */
      const g = t * (nres - 1), i0 = Math.min(nres - 1, Math.floor(g));
      const i1 = Math.min(nres - 1, i0 + 1), fr = smoothstep(g - i0);
      const k0 = ss[a + i0], k1 = ss[a + i1];
      let hw = HW[k0] + (HW[k1] - HW[k0]) * fr;
      let ht = HT[k0] + (HT[k1] - HT[k0]) * fr;
      /* taper the open ends so a run finishes in a point, not a stub */
      const edge = Math.min(s, steps - s) / (SUB * 0.8);
      const taper = 0.22 + 0.78 * smoothstep(Math.min(1, edge));
      hw *= taper; ht *= taper;

      const alpha = alphaAt(P.x, P.y, P.z);
      for (let k = 0; k < PROFILE; k++) {
        const du = SEC.uu[k] * hw, dvv = SEC.vv[k] * ht;
        pos.push(P.x + W.x * du + N.x * dvv,
                 P.y + W.y * du + N.y * dvv,
                 P.z + W.z * du + N.z * dvv);
        /* outward normal: the implicit gradient, scaled back out of the
           anisotropic cross-section so a flat ribbon still lights correctly */
        tmp.copy(W).multiplyScalar(SEC.nu[k] / hw).addScaledVector(N, SEC.nv[k] / ht);
        if (tmp.lengthSq() < 1e-12) tmp.copy(N);
        tmp.normalize();
        nrm.push(tmp.x, tmp.y, tmp.z);
        col.push(1, 1, 1, alpha);
      }

      /* Decide the winding once per run by checking the first quad against
         the normal we already know, instead of assuming a handedness. */
      if (s === 1 && !haveWinding) {
        const i00 = ring0 * 3, i01 = (ring0 + 1) * 3, i10 = (ring0 + PROFILE) * 3;
        triA.set(pos[i01] - pos[i00], pos[i01 + 1] - pos[i00 + 1], pos[i01 + 2] - pos[i00 + 2]);
        triB.set(pos[i10] - pos[i00], pos[i10 + 1] - pos[i00 + 1], pos[i10 + 2] - pos[i00 + 2]);
        triN.crossVectors(triA, triB);
        flip = triN.dot(new THREE.Vector3(nrm[i00], nrm[i00 + 1], nrm[i00 + 2])) < 0;
        haveWinding = true;
      }
    }

    const quad = (p0, p1, p2, p3) => {
      if (flip) idx.push(p0, p2, p1, p0, p3, p2);
      else idx.push(p0, p1, p2, p0, p2, p3);
    };
    for (let s = 0; s < steps; s++) {
      const r0 = ring0 + s * PROFILE, r1 = r0 + PROFILE;
      for (let k = 0; k < PROFILE; k++) {
        const k1 = (k + 1) % PROFILE;
        quad(r0 + k, r0 + k1, r1 + k1, r1 + k);
      }
    }

    /* caps: a fan around the centre of each end ring, flat along the path */
    for (const [ringBase, sign] of [[ring0, -1], [ring0 + steps * PROFILE, 1]]) {
      const cIdx = pos.length / 3;
      let cx = 0, cy = 0, cz = 0, ca4 = 0;
      for (let k = 0; k < PROFILE; k++) {
        const j = (ringBase + k) * 3;
        cx += pos[j] / PROFILE; cy += pos[j + 1] / PROFILE; cz += pos[j + 2] / PROFILE;
        ca4 += col[(ringBase + k) * 4 + 3] / PROFILE;
      }
      cA.getTangent(sign < 0 ? 0 : 1, T).normalize().multiplyScalar(sign);
      pos.push(cx, cy, cz); nrm.push(T.x, T.y, T.z); col.push(1, 1, 1, ca4);
      /* the rim is duplicated so the cap keeps its own flat normal */
      const rim0 = pos.length / 3;
      for (let k = 0; k < PROFILE; k++) {
        const j = (ringBase + k) * 3;
        pos.push(pos[j], pos[j + 1], pos[j + 2]);
        nrm.push(T.x, T.y, T.z);
        col.push(1, 1, 1, col[(ringBase + k) * 4 + 3]);
      }
      for (let k = 0; k < PROFILE; k++) {
        const k1 = (k + 1) % PROFILE;
        if (sign < 0) idx.push(cIdx, rim0 + k1, rim0 + k);
        else idx.push(cIdx, rim0 + k, rim0 + k1);
      }
    }
  }

  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
  g.setAttribute('normal', new THREE.Float32BufferAttribute(nrm, 3));
  /* 4 components, so the alpha reaches the material -- this is USE_COLOR_ALPHA
     and needs vertexColors and transparent on a plain (non-instanced) mesh */
  g.setAttribute('color', new THREE.Float32BufferAttribute(col, 4));
  g.setIndex(idx);
  g.computeBoundingSphere();
  return g;
}
