/* shepherd-score alignment viewer --------------------------------------------
   Plays back a recorded single-start alignment: a reference molecule stays put
   while a fit molecule is moved rigidly by the optimizer, one Adam step per
   frame, maximising one of three overlap scores.

     vol_noH  Gaussian volume overlap of heavy atoms     (atoms only, no H)
     surf     Gaussian overlap of surface point clouds   (atoms + surface points)
     pharm    vector-aware pharmacophore overlap         (atoms + pharmacophores)

   The fit is stored once, in its start pose, and each frame is the rigid
   transform from that pose (quaternion xyzw + translation), so playback is a
   matrix update rather than new geometry.  Between recorded steps the pose is
   slerped, which keeps slow playback smooth.

   Data comes from tools/build_alignment.py.
   Mount with:  mount(boxEl, window.SHEPHERD_ALIGN, { mode: 'vol_noH' })
---------------------------------------------------------------------------- */
import * as THREE from 'three';
/* Trackball rather than Orbit, as in traj-viewer.js: free tumbling, no pole lock. */
import { TrackballControls } from 'three/addons/controls/TrackballControls.js';

/* ---------- styling ---------- */
const ELEM = {
  H:  { c: 0xe8e8ea, r: 0.30 }, C:  { c: 0x4c5057, r: 0.46 },
  N:  { c: 0x3b6fd4, r: 0.44 }, O:  { c: 0xd4453b, r: 0.43 },
  F:  { c: 0x5bbf7a, r: 0.40 }, P:  { c: 0xe08a2e, r: 0.52 },
  S:  { c: 0xd9b63a, r: 0.52 }, Cl: { c: 0x4fae5a, r: 0.50 },
  Br: { c: 0xa2543a, r: 0.55 }, I:  { c: 0x8c4fc4, r: 0.60 },
};
const FALLBACK = { c: 0x9aa0a8, r: 0.45 };
/* The two molecules differ only in their carbons: the reference is a pale,
   quiet grey and the moving fit is violet, a hue no element or pharmacophore
   on the page uses. */
export const REF_C = 0xa9aeb6;
export const FIT_C = 0x8a6fd8;
export const SURF_REF = 0xa9aeb6;
export const SURF_FIT = 0x8a6fd8;

/* shepherd_score pharmacophore palette, indexed like P_TYPES.  Hydrophobe is
   lifted from (0.2,0.2,0.2) so a translucent sphere survives a dark page. */
export const PHARM = [
  { name: 'Acceptor',   c: [0.62, 0.03, 0.35] },
  { name: 'Donor',      c: [0.00, 0.55, 0.55] },
  { name: 'Aromatic',   c: [0.85, 0.50, 0.00] },
  { name: 'Hydrophobe', c: [0.42, 0.44, 0.47], r: 0.5 },
  { name: 'Halogen',    c: [0.00, 0.80, 0.10] },
  { name: 'Cation',     c: [0.00, 0.00, 1.00], r: 0.5 },
  { name: 'Anion',      c: [1.00, 0.00, 0.00], r: 0.5 },
  { name: 'ZnBinder',   c: [1.00, 0.50, 0.50], r: 0.5 },
];
const VECTOR_TYPES = new Set([0, 1, 2, 4]);   /* drawn as arrows, not spheres */
const VEC_SCALE = 2;                           /* as in the ChimeraX export */

const MODES = ['vol_noH', 'surf', 'pharm'];

export function mount(box, data, opts = {}) {
  const stage = box.querySelector('.tv-stage');
  stage.querySelectorAll('.tv-fallback').forEach(n => n.remove());

  const vec = (a) => new THREE.Vector3(a[0], a[1], a[2]);
  const poses = {};
  for (const m of MODES) {
    const d = data.modes[m];
    poses[m] = {
      q: d.quat.map(q => new THREE.Quaternion(q[0], q[1], q[2], q[3])),
      t: d.trans.map(vec),
      s: d.scores,
      n: d.scores.length,
    };
  }

  /* ---------- framing ----------
     Bounding sphere over everything that is ever on screen: the reference with
     its surface, and the fit at its start pose and at each mode's end pose. */
  let fitR = 1;
  const grow = (p) => { fitR = Math.max(fitR, p.length()); };
  for (const p of [...data.ref.pos, ...data.ref.surf]) grow(vec(p));
  for (const m of MODES) {
    const P = poses[m], last = P.n - 1;
    for (const f of [0, last]) for (const p of [...data.fit.pos, ...data.fit.surf]) {
      grow(vec(p).applyQuaternion(P.q[f]).add(P.t[f]));
    }
  }
  const FOV = 38;
  const DIST = (fitR * 0.86) / Math.tan((FOV / 2) * Math.PI / 180);

  /* ---------- opening viewpoint ----------
     Look from the direction that spreads the aligned heavy atoms furthest
     apart on screen (least crowding), then roll so the long axis of the
     projection lies across the frame. */
  const { HOME, HOME_UP } = (() => {
    const pts = [];
    data.ref.pos.forEach((p, i) => { if (data.ref.elements[i] !== 'H') pts.push(vec(p)); });
    const V = poses.vol_noH, last = V.n - 1;
    data.fit.pos.forEach((p, i) => {
      if (data.fit.elements[i] !== 'H') pts.push(vec(p).applyQuaternion(V.q[last]).add(V.t[last]));
    });
    let best = null, bestC = Infinity;
    const golden = Math.PI * (3 - Math.sqrt(5)), N = 400, a = new THREE.Vector3();
    for (let k = 0; k < N; k++) {           /* Fibonacci hemisphere */
      const z = 1 - k / N, rr = Math.sqrt(1 - z * z), th = k * golden;
      const dir = new THREE.Vector3(Math.cos(th) * rr, Math.sin(th) * rr, z);
      let c = 0;
      for (let i = 0; i < pts.length; i++) for (let j = i + 1; j < pts.length; j++) {
        a.subVectors(pts[i], pts[j]);
        const d2 = a.lengthSq() - a.dot(dir) ** 2;
        c += Math.exp(-d2 / 0.5);
      }
      if (c < bestC) { bestC = c; best = dir; }
    }
    /* long axis of the projection, by 2D PCA in the view plane */
    const u = new THREE.Vector3().crossVectors(best, Math.abs(best.y) < 0.9
      ? new THREE.Vector3(0, 1, 0) : new THREE.Vector3(1, 0, 0)).normalize();
    const w = new THREE.Vector3().crossVectors(best, u);
    let sxx = 0, syy = 0, sxy = 0;
    for (const p of pts) { const x = p.dot(u), y = p.dot(w); sxx += x * x; syy += y * y; sxy += x * y; }
    const ang = 0.5 * Math.atan2(2 * sxy, sxx - syy);
    const longAxis = u.clone().multiplyScalar(Math.cos(ang)).addScaledVector(w, Math.sin(ang));
    return { HOME: best.clone().multiplyScalar(DIST),
             HOME_UP: new THREE.Vector3().crossVectors(best, longAxis).normalize() };
  })();

  /* ---------- scene ---------- */
  const scene = new THREE.Scene();
  const camera = new THREE.PerspectiveCamera(FOV, 1, 0.5, 2000);
  camera.position.copy(HOME); camera.up.copy(HOME_UP); camera.lookAt(0, 0, 0);

  const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
  renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
  renderer.outputColorSpace = THREE.SRGBColorSpace;
  stage.appendChild(renderer.domElement);
  renderer.domElement.classList.add('tv-canvas');

  /* camera-attached light rig (see traj-viewer.js) */
  scene.add(camera);
  const rig = new THREE.Group();
  camera.add(rig);
  const addLight = (color, intensity, pos) => {
    const l = new THREE.DirectionalLight(color, intensity);
    l.position.set(...pos);
    l.target.position.set(0, 0, -1);
    rig.add(l, l.target);
  };
  addLight(0xffffff, 1.55, [0.5, 0.7, 0.3]);
  addLight(0xcfdcf2, 0.75, [-0.85, -0.4, 0.25]);
  addLight(0xffffff, 0.5, [0, 0.35, -2.5]);
  scene.add(new THREE.HemisphereLight(0xf4f6f9, 0x9aa1aa, 1.15));

  const mkMat = (o) => new THREE.MeshStandardMaterial(
    Object.assign({ roughness: 0.42, metalness: 0.03 }, o));
  const disposables = [];
  const inst = (geo, mat, n) => {
    const m = new THREE.InstancedMesh(geo, mat, Math.max(1, n));
    m.count = n;
    m.frustumCulled = false;
    disposables.push(m);
    return m;
  };

  const M = new THREE.Matrix4(), Q = new THREE.Quaternion(), C = new THREE.Color();
  const UPV = new THREE.Vector3(0, 1, 0), ONE = new THREE.Vector3(1, 1, 1);
  function cylinder(mesh, i, a, b, r) {
    const d = new THREE.Vector3().subVectors(b, a);
    const len = d.length() || 1e-6;
    Q.setFromUnitVectors(UPV, d.clone().divideScalar(len));
    M.compose(a.clone().addScaledVector(d, 0.5), Q, new THREE.Vector3(r, len, r));
    mesh.setMatrixAt(i, M);
  }

  /* One molecule -> { heavy, hydrogens, surf, pharm } groups, built once. */
  function buildMolecule(mol, carbon, surfColor, isFit) {
    const P = mol.pos.map(vec);
    const style = mol.elements.map(e => (e === 'C' ? { c: carbon, r: ELEM.C.r } : ELEM[e] || FALLBACK));
    const isH = mol.elements.map(e => e === 'H');

    const out = { heavy: new THREE.Group(), hydrogens: new THREE.Group(),
                  surf: new THREE.Group(), pharm: new THREE.Group() };

    for (const wantH of [false, true]) {
      const idx = mol.elements.map((_, i) => i).filter(i => isH[i] === wantH);
      const atoms = inst(new THREE.SphereGeometry(1, 22, 16), mkMat({}), idx.length);
      idx.forEach((i, k) => {
        const s = style[i].r * 0.55;
        M.compose(P[i], Q.identity(), new THREE.Vector3(s, s, s));
        atoms.setMatrixAt(k, M);
        atoms.setColorAt(k, C.setHex(style[i].c));
      });
      /* half-bonds, each coloured by its own atom; a bond touching H goes with H */
      const bl = mol.bonds.filter(([a, b]) => (isH[a] || isH[b]) === wantH);
      const bonds = inst(new THREE.CylinderGeometry(0.11, 0.11, 1, 12), mkMat({}), bl.length * 2);
      bl.forEach(([a, b], k) => {
        const mid = P[a].clone().add(P[b]).multiplyScalar(0.5);
        cylinder(bonds, 2 * k, P[a], mid, 1);
        cylinder(bonds, 2 * k + 1, mid, P[b], 1);
        bonds.setColorAt(2 * k, C.setHex(style[a].c));
        bonds.setColorAt(2 * k + 1, C.setHex(style[b].c));
      });
      (wantH ? out.hydrogens : out.heavy).add(atoms, bonds);
    }

    const sp = inst(new THREE.SphereGeometry(1, 10, 8),
      mkMat({ color: surfColor, roughness: 0.6, transparent: true, opacity: isFit ? 0.9 : 0.6 }),
      mol.surf.length);
    mol.surf.forEach((p, k) => { M.compose(vec(p), Q.identity(), new THREE.Vector3(0.11, 0.11, 0.11)); sp.setMatrixAt(k, M); });
    out.surf.add(sp);

    /* pharmacophores: vectors as arrows of length 2 (aromatic both ways),
       the rest as translucent spheres.  The fit's are opaque, the reference's
       faded, so the two sets can be told apart where they meet. */
    const { types, ancs, vecs } = mol.pharm;
    const arrows = [], spheres = [];
    types.forEach((t, i) => {
      const a = vec(ancs[i]), v = vec(vecs[i]);
      if (VECTOR_TYPES.has(t) && v.lengthSq() > 1e-8) {
        arrows.push([t, a, v]);
        if (t === 2) arrows.push([t, a, v.clone().negate()]);
      } else spheres.push([t, a]);
    });
    const op = isFit ? 0.95 : 0.5;
    const shafts = inst(new THREE.CylinderGeometry(0.1, 0.1, 1, 10),
      mkMat({ transparent: !isFit, opacity: op }), arrows.length);
    const heads = inst(new THREE.ConeGeometry(0.22, 0.48, 14),
      mkMat({ transparent: !isFit, opacity: op }), arrows.length);
    arrows.forEach(([t, a, v], k) => {
      const tip = a.clone().addScaledVector(v, VEC_SCALE);
      const base = tip.clone().addScaledVector(v, -0.48);
      cylinder(shafts, k, a, base, 1);
      Q.setFromUnitVectors(UPV, v.clone().normalize());
      M.compose(base.clone().addScaledVector(v, 0.24), Q, ONE);
      heads.setMatrixAt(k, M);
      C.setRGB(...PHARM[t].c, THREE.SRGBColorSpace);
      shafts.setColorAt(k, C); heads.setColorAt(k, C);
    });
    const balls = inst(new THREE.SphereGeometry(1, 22, 16),
      mkMat({ transparent: true, opacity: isFit ? 0.55 : 0.3, depthWrite: false }), spheres.length);
    spheres.forEach(([t, a], k) => {
      const r = PHARM[t]?.r ?? 0.5;
      M.compose(a, Q.identity(), new THREE.Vector3(r, r, r));
      balls.setMatrixAt(k, M);
      balls.setColorAt(k, C.setRGB(...(PHARM[t]?.c ?? [0.5, 0.5, 0.5]), THREE.SRGBColorSpace));
    });
    out.pharm.add(shafts, heads, balls);

    const g = new THREE.Group();
    g.add(out.heavy, out.hydrogens, out.surf, out.pharm);
    out.group = g;
    return out;
  }

  const root = new THREE.Group();
  scene.add(root);
  const ref = buildMolecule(data.ref, REF_C, SURF_REF, false);
  const fit = buildMolecule(data.fit, FIT_C, SURF_FIT, true);
  root.add(ref.group, fit.group);

  const controls = new TrackballControls(camera, renderer.domElement);
  controls.rotateSpeed = 3.2;
  controls.zoomSpeed = 1.1;
  controls.noPan = true;
  controls.dynamicDampingFactor = 0.14;
  controls.minDistance = DIST * 0.3;
  controls.maxDistance = DIST * 3;
  controls.keys = [];

  /* ---------- UI ---------- */
  const ui = document.createElement('div');
  ui.className = 'tv-ui';
  ui.innerHTML = `
    <div class="tv-row tv-transport">
      <button class="tv-btn tv-play" type="button" aria-label="Play">
        <span class="tv-ico-play" aria-hidden="true"></span><span class="tv-play-label">Play</span>
      </button>
      <input class="tv-slider" type="range" min="0" value="0" step="1" aria-label="Optimizer step">
      <span class="tv-step" aria-live="off"></span>
    </div>
    <div class="av-score">
      <div class="av-score-head">
        <span class="av-score-label">Similarity</span>
        <span class="av-score-val">0.000</span>
        <span class="av-score-best"></span>
      </div>
      <div class="av-spark" role="img">
        <svg viewBox="0 0 1000 100" preserveAspectRatio="none" aria-hidden="true">
          <line class="av-grid" x1="0" x2="1000" y1="100" y2="100"/>
          <line class="av-best" x1="0" x2="1000"/>
          <polyline class="av-line" fill="none"/>
          <line class="av-cursor" y1="0" y2="100"/>
        </svg>
        <span class="av-dot"></span>
        <span class="av-tip" hidden></span>
      </div>
    </div>
    <div class="tv-row tv-toggles">
      <select class="tv-speed" aria-label="Playback speed">
        <option value="0.5">0.5&times;</option>
        <option value="1" selected>1&times;</option>
        <option value="2">2&times;</option>
      </select>
      <button class="tv-btn tv-reset" type="button">Reset view</button>
    </div>`;
  box.appendChild(ui);
  const $ = (s) => ui.querySelector(s);
  const slider = $('.tv-slider'), playBtn = $('.tv-play'), stepOut = $('.tv-step'),
        speedSel = $('.tv-speed'), scoreVal = $('.av-score-val'), bestOut = $('.av-score-best'),
        spark = $('.av-spark'), line = $('.av-line'), bestLine = $('.av-best'),
        cursor = $('.av-cursor'), dot = $('.av-dot'), tip = $('.av-tip');

  /* ---------- state ---------- */
  let mode = MODES.includes(opts.mode) ? opts.mode : 'vol_noH';
  let pos = 0;              /* fractional frame */
  let playing = false, speed = 1, last = 0;
  const FPS = 30;
  const P = () => poses[mode];
  const yOf = (s) => 100 - s * 100;   /* score axis is fixed 0..1 for every mode */

  function applyLayers() {
    for (const m of [ref, fit]) {
      m.hydrogens.visible = mode !== 'vol_noH';
      m.surf.visible = mode === 'surf';
      m.pharm.visible = mode === 'pharm';
    }
  }

  function drawSpark() {
    const { s, n } = P();
    line.setAttribute('points', s.map((v, i) => `${(i / (n - 1)) * 1000},${yOf(v)}`).join(' '));
    const best = data.info.library_best[mode];
    bestLine.setAttribute('y1', yOf(best)); bestLine.setAttribute('y2', yOf(best));
    bestOut.textContent = `best of 50 starts ${best.toFixed(3)}`;
    spark.setAttribute('aria-label',
      `Score by optimizer step: ${s[0].toFixed(2)} at the start, ${s[n - 1].toFixed(2)} after ${n - 1} steps`);
  }

  const qTmp = new THREE.Quaternion();
  function setPos(p, fromSlider) {
    const { q, t, s, n } = P();
    pos = Math.max(0, Math.min(n - 1, p));
    const i = Math.min(n - 2, Math.floor(pos)), f = pos - i;
    qTmp.slerpQuaternions(q[i], q[i + 1], f);
    fit.group.quaternion.copy(qTmp);
    fit.group.position.lerpVectors(t[i], t[i + 1], f);
    const sc = s[i] + (s[i + 1] - s[i]) * f;
    const k = Math.round(pos);
    if (!fromSlider) slider.value = k;
    stepOut.textContent = `step ${k} / ${n - 1}`;
    scoreVal.textContent = sc.toFixed(3);
    const x = (pos / (n - 1)) * 100;
    cursor.setAttribute('x1', x * 10); cursor.setAttribute('x2', x * 10);
    dot.style.left = `${x}%`;
    dot.style.top = `${yOf(sc)}%`;
  }

  function setPlaying(p) {
    playing = p;
    playBtn.classList.toggle('is-playing', p);
    playBtn.setAttribute('aria-label', p ? 'Pause' : 'Play');
    $('.tv-play-label').textContent = p ? 'Pause' : 'Play';
  }

  const reduced = matchMedia('(prefers-reduced-motion: reduce)').matches;
  function setMode(m, { play = !reduced } = {}) {
    if (!MODES.includes(m)) return;
    mode = m;
    slider.max = P().n - 1;
    applyLayers();
    drawSpark();
    setPos(0);
    setPlaying(play);
  }

  playBtn.addEventListener('click', () => {
    if (!playing && pos >= P().n - 1) setPos(0);
    setPlaying(!playing);
  });
  slider.addEventListener('input', () => { setPlaying(false); setPos(+slider.value, true); });
  speedSel.addEventListener('change', () => { speed = +speedSel.value; });
  $('.tv-reset').addEventListener('click', () => {
    controls.reset();
    camera.position.copy(HOME); camera.up.copy(HOME_UP); camera.lookAt(0, 0, 0);
  });

  /* The score chart doubles as a scrubber: hover reads a step off it, press
     and drag seeks there. */
  const stepAt = (e) => {
    const r = spark.getBoundingClientRect();
    return Math.round(Math.max(0, Math.min(1, (e.clientX - r.left) / r.width)) * (P().n - 1));
  };
  const showTip = (k) => {
    const { s, n } = P();
    tip.hidden = false;
    tip.textContent = `step ${k} · ${s[k].toFixed(3)}`;
    const x = (k / (n - 1)) * 100;
    tip.style.left = `${x}%`;
    tip.style.transform = `translateX(${x > 70 ? '-100%' : x < 30 ? '0' : '-50%'})`;
  };
  let dragging = false;
  spark.addEventListener('pointermove', (e) => {
    const k = stepAt(e);
    showTip(k);
    if (dragging) setPos(k);
  });
  spark.addEventListener('pointerdown', (e) => {
    dragging = true; spark.setPointerCapture(e.pointerId);
    setPlaying(false); setPos(stepAt(e));
  });
  spark.addEventListener('pointerup', () => { dragging = false; });
  spark.addEventListener('pointerleave', () => { if (!dragging) tip.hidden = true; });

  box.addEventListener('keydown', (e) => {
    if (e.target.closest('select, input, button')) return;
    if (e.key === ' ') { e.preventDefault(); playBtn.click(); }
    else if (e.key === 'ArrowRight') { setPlaying(false); setPos(Math.round(pos) + 1); }
    else if (e.key === 'ArrowLeft') { setPlaying(false); setPos(Math.round(pos) - 1); }
  });
  box.tabIndex = 0;

  /* ---------- resize + loop ---------- */
  function resize() {
    const w = stage.clientWidth, h = stage.clientHeight;
    if (!w || !h) return;
    renderer.setSize(w, h, false);
    camera.aspect = w / h;
    camera.updateProjectionMatrix();
    controls.handleResize();
  }
  const ro = new ResizeObserver(resize); ro.observe(stage);
  resize();

  let visible = true;
  const vo = new IntersectionObserver((es) => { visible = es[0].isIntersecting; }, { threshold: 0.01 });
  vo.observe(stage);

  let alive = true;
  function tick(now) {
    if (!alive) return;
    requestAnimationFrame(tick);
    const dt = Math.min(0.1, (now - last) / 1000) || 0; last = now;
    if (!visible) return;
    if (playing) {
      const end = P().n - 1;
      const next = pos + dt * FPS * speed;
      if (next >= end) { setPos(end); setPlaying(false); } else setPos(next);
    }
    controls.update();
    renderer.render(scene, camera);
  }

  setMode(mode, { play: false });
  requestAnimationFrame((t) => { last = t; tick(t); });

  /* start playing the first time the viewer is mostly on screen */
  if (!reduced && opts.autoplay !== false) {
    const io = new IntersectionObserver((es) => {
      if (es[0].isIntersecting) { setPlaying(true); io.disconnect(); }
    }, { threshold: 0.4 });
    io.observe(stage);
  }

  function dispose() {
    alive = false;
    ro.disconnect(); vo.disconnect();
    controls.dispose();
    for (const m of disposables) { m.geometry.dispose(); m.material.dispose(); m.dispose(); }
    renderer.dispose();
    renderer.domElement.remove();
    ui.remove();
  }

  const seek = (k) => { setPlaying(false); setPos(k); };
  return { setMode, setPlaying, seek, dispose, get mode() { return mode; } };
}
