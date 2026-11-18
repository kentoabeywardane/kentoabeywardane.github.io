/* Choosing an opening camera ---------------------------------------------------
   A molecule comes out of the dump in an arbitrary orientation, and looking at
   a flat molecule edge-on hides everything that makes it readable.  These are
   the pieces the viewers use to pick a first view from the geometry itself,
   rather than from hand-tuned angles that would not survive a new trajectory.

   Measured on the profile trajectory, the closest pair of projected atoms sits
   0.22 A apart looking straight down the best-fit plane normal and 0.47 A at a
   22 degree tilt -- spreading the atoms out and tilting for a three-quarter
   view want the same thing, they do not trade off against each other.
---------------------------------------------------------------------------- */
import * as THREE from 'three';

/* Principal axes of a point set, largest variance first, plus the points
   recentred on their own mean.  Jacobi rotation: the matrix is a symmetric
   3x3, so this converges in a handful of sweeps. */
export function principalAxes(pts) {
  if (pts.length < 3) return null;
  const m = [0, 0, 0];
  for (const p of pts) for (let k = 0; k < 3; k++) m[k] += p[k] / pts.length;
  const c = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (const p of pts)
    for (let a = 0; a < 3; a++)
      for (let b = 0; b < 3; b++) c[a][b] += (p[a] - m[a]) * (p[b] - m[b]) / pts.length;

  const v = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  for (let sweep = 0; sweep < 24; sweep++) {
    let p = 0, q = 1, off = Math.abs(c[0][1]);
    if (Math.abs(c[0][2]) > off) { p = 0; q = 2; off = Math.abs(c[0][2]); }
    if (Math.abs(c[1][2]) > off) { p = 1; q = 2; off = Math.abs(c[1][2]); }
    if (off < 1e-10) break;
    const theta = (c[q][q] - c[p][p]) / (2 * c[p][q]);
    const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
    const cs = 1 / Math.sqrt(t * t + 1), sn = t * cs;
    for (let k = 0; k < 3; k++) {
      const ckp = c[k][p], ckq = c[k][q];
      c[k][p] = cs * ckp - sn * ckq; c[k][q] = sn * ckp + cs * ckq;
    }
    for (let k = 0; k < 3; k++) {
      const cpk = c[p][k], cqk = c[q][k];
      c[p][k] = cs * cpk - sn * cqk; c[q][k] = sn * cpk + cs * cqk;
      const vkp = v[k][p], vkq = v[k][q];
      v[k][p] = cs * vkp - sn * vkq; v[k][q] = sn * vkp + cs * vkq;
    }
  }
  const order = [0, 1, 2].sort((a, b) => c[b][b] - c[a][a]);
  const axis = (i) => new THREE.Vector3(v[0][order[i]], v[1][order[i]], v[2][order[i]]).normalize();
  return {
    major: axis(0), mid: axis(1), normal: axis(2),
    pts: pts.map((p) => new THREE.Vector3(p[0] - m[0], p[1] - m[1], p[2] - m[2])),
    center: new THREE.Vector3(m[0], m[1], m[2]),
  };
}

/* How badly a view direction piles the atoms on top of each other.  Sum of
   1/(projected distance) over every pair: it blows up when two atoms land in
   the same place, so the minimum is the least-occluded direction. */
export function crowding(pts, dir) {
  const e1 = new THREE.Vector3(), e2 = new THREE.Vector3();
  e1.set(0, 0, 1);
  if (Math.abs(dir.z) > 0.9) e1.set(1, 0, 0);
  e1.cross(dir).normalize();
  e2.crossVectors(dir, e1).normalize();
  const n = pts.length, a = new Float64Array(n), b = new Float64Array(n);
  for (let i = 0; i < n; i++) { a[i] = pts[i].dot(e1); b[i] = pts[i].dot(e2); }
  let sum = 0;
  for (let i = 0; i < n; i++)
    for (let j = i + 1; j < n; j++) {
      const dx = a[i] - a[j], dy = b[i] - b[j];
      sum += 1 / (Math.sqrt(dx * dx + dy * dy) + 0.25);
    }
  return sum;
}

/* Least-squares fit of perp = c0 + c1*t + c2*t^2.  Only the sign of c2 is
   used: it separates a molecule that sags in the middle from one that arches. */
export function curvature(t, perp) {
  const M = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]];
  for (let i = 0; i < t.length; i++) {
    const p = [1, t[i], t[i] * t[i]];
    for (let r = 0; r < 3; r++) {
      for (let c = 0; c < 3; c++) M[r][c] += p[r] * p[c];
      M[r][3] += p[r] * perp[i];
    }
  }
  for (let col = 0; col < 3; col++) {              /* Gaussian elimination */
    let piv = col;
    for (let r = col + 1; r < 3; r++) if (Math.abs(M[r][col]) > Math.abs(M[piv][col])) piv = r;
    if (Math.abs(M[piv][col]) < 1e-12) return 0;
    [M[col], M[piv]] = [M[piv], M[col]];
    for (let r = 0; r < 3; r++) {
      if (r === col) continue;
      const k = M[r][col] / M[col][col];
      for (let c = col; c < 4; c++) M[r][c] -= k * M[col][c];
    }
  }
  return M[2][3] / M[2][2];
}

/* The least-crowded direction within `maxTilt` of `anchor`, searched on an
   evenly sampled spherical cap (polar angle from the cosine so the samples do
   not bunch at the pole, azimuth by the golden angle). */
export function leastCrowded(pts, anchor, maxTilt, axes, samples = 700) {
  const dir = anchor.clone().normalize();
  let best = crowding(pts, dir);
  const cand = new THREE.Vector3();
  const GOLD = Math.PI * (1 + Math.sqrt(5));
  const u = axes?.major ?? new THREE.Vector3(1, 0, 0);
  const w = axes?.mid ?? new THREE.Vector3(0, 1, 0);
  /* a basis square to the anchor, whatever was passed in */
  const e1 = u.clone().addScaledVector(dir, -u.dot(dir));
  if (e1.lengthSq() < 1e-8) e1.set(dir.y, -dir.x, 0);
  e1.normalize();
  const e2 = new THREE.Vector3().crossVectors(dir, e1).normalize();
  void w;
  for (let i = 0; i < samples; i++) {
    const cz = 1 - (i / samples) * (1 - Math.cos(maxTilt));
    const sr = Math.sqrt(Math.max(0, 1 - cz * cz)), az = GOLD * i;
    cand.copy(dir).multiplyScalar(cz)
        .addScaledVector(e1, sr * Math.cos(az))
        .addScaledVector(e2, sr * Math.sin(az));
    cand.normalize();
    const s = crowding(pts, cand);
    if (s < best) { best = s; dir.copy(cand); }
  }
  return dir;
}

/* The roll that lays the long axis across the frame and leaves the molecule
   sagging rather than arching -- ends high, middle low.  Returns the camera
   up-vector for the given view direction. */
export function uprightFor(pts, dir) {
  const e1 = new THREE.Vector3(0, 0, 1);
  if (Math.abs(dir.z) > 0.9) e1.set(1, 0, 0);
  e1.cross(dir).normalize();
  const e2 = new THREE.Vector3().crossVectors(dir, e1).normalize();
  let sxx = 0, syy = 0, sxy = 0;
  for (const p of pts) {
    const a = p.dot(e1), b = p.dot(e2);
    sxx += a * a; syy += b * b; sxy += a * b;
  }
  const th = 0.5 * Math.atan2(2 * sxy, sxx - syy);      /* 2x2 PCA, closed form */
  const right = e1.clone().multiplyScalar(Math.cos(th))
                  .addScaledVector(e2, Math.sin(th)).normalize();
  const up = new THREE.Vector3().crossVectors(dir, right).normalize();
  const t = [], perp = [];
  for (const p of pts) { t.push(p.dot(right)); perp.push(p.dot(up)); }
  if (curvature(t, perp) < 0) up.negate();              /* arching -> flip */
  return up;
}
