/* ShEPhERD trajectory viewer -------------------------------------------------
   Renders a precomputed reverse-diffusion trajectory: the molecule (x1) and the
   interaction profile it is generated against -- the ESP surface (x3) and the
   pharmacophores (x4).  Data comes from tools/build_traj.py.

   Mount with:  TrajViewer.mount(containerEl, window.SHEPHERD2_TRAJ)
---------------------------------------------------------------------------- */
import * as THREE from 'three';
/* Trackball rather than Orbit: OrbitControls pins an up-axis, so the view
   locks up at the poles and the molecule cannot be tumbled freely. */
import { TrackballControls } from 'three/addons/controls/TrackballControls.js';

/* ---------- element styling (CPK, slightly desaturated) ---------- */
const ELEM = {
  H:  { c: 0xe8e8ea, r: 0.30 }, C:  { c: 0x4c5057, r: 0.46 },
  N:  { c: 0x3b6fd4, r: 0.44 }, O:  { c: 0xd4453b, r: 0.43 },
  F:  { c: 0x5bbf7a, r: 0.40 }, P:  { c: 0xe08a2e, r: 0.52 },
  S:  { c: 0xd9b63a, r: 0.52 }, Cl: { c: 0x4fae5a, r: 0.50 },
  Br: { c: 0xa2543a, r: 0.55 }, I:  { c: 0x8c4fc4, r: 0.60 },
  Si: { c: 0xb0a090, r: 0.56 },
};
const FALLBACK = { c: 0x9aa0a8, r: 0.45 };
/* covalent radii (A) for the distance-based bond test */
const COV = { H: 0.31, C: 0.76, N: 0.71, O: 0.66, F: 0.57, P: 1.07, S: 1.05,
              Cl: 1.02, Br: 1.20, I: 1.39, Si: 1.11 };

/* ---------- helpers ---------- */
function unb64(s, Type) {
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
function espColor(t, out) {
  const a = Math.abs(t);
  const lo = [0.88, 0.88, 0.86];
  const hi = t >= 0 ? [0.11, 0.38, 0.85] : [0.84, 0.17, 0.14];
  const k = a;
  out.setRGB(lo[0] + (hi[0] - lo[0]) * k,
             lo[1] + (hi[1] - lo[1]) * k,
             lo[2] + (hi[2] - lo[2]) * k, THREE.SRGBColorSpace);
  return out;
}

export function mount(box, data, opts = {}) {
  const stage = box.querySelector('.tv-stage') || box.appendChild(
    Object.assign(document.createElement('div'), { className: 'tv-stage' }));
  stage.querySelectorAll('.tv-fallback').forEach(n => n.remove());
  const host = stage;
  const S = data.posScale;
  const nF = data.nframes;

  /* ---------- decode ---------- */
  const x1c = unb64(data.x1.counts, Uint16Array);
  const x1o = offsets(x1c);
  const x1e = unb64(data.x1.elem, Uint8Array);
  const x1p = unb64(data.x1.pos, Int16Array);

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

  const palette = data.colors.map(
    c => new THREE.Color().setRGB(c[0], c[1], c[2], THREE.SRGBColorSpace));
  const elemStyle = data.elements.map(e => ELEM[e] || FALLBACK);
  const elemCov = data.elements.map(e => COV[e] ?? 0.9);

  const MAXA = Math.max(...x1c);
  const MAXS = Math.max(...s4c);
  const MAXR = Math.max(...a4c);
  const MAXB = MAXA * 4;

  /* ---------- scene ---------- */
  /* frame the camera on the settled structure of the last frame */
  let fitR = 1;
  {
    const f = nF - 1, a0 = x1o[f];
    for (let i = 0; i < x1c[f]; i++)
      for (let k = 0; k < 3; k++) fitR = Math.max(fitR, Math.abs(x1p[(a0 + i) * 3 + k] / S));
    for (let i = 0; i < n3; i++)
      for (let k = 0; k < 3; k++) fitR = Math.max(fitR, Math.abs(x3p[(f * n3 + i) * 3 + k] / S));
  }
  const FOV = 38;
  const DIST = (fitR * 1.45) / Math.tan((FOV / 2) * Math.PI / 180);

  const scene = new THREE.Scene();
  const camera = new THREE.PerspectiveCamera(FOV, 1, 0.5, 2000);
  camera.position.set(0, 0, DIST);

  const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
  renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
  renderer.outputColorSpace = THREE.SRGBColorSpace;
  host.appendChild(renderer.domElement);
  renderer.domElement.classList.add('tv-canvas');

  /* Lighting rides on the camera, so whichever way the molecule is turned the
     face you are looking at is the face that is lit.  World-fixed lights leave
     a dark side as soon as you rotate past them.  Directional lights aim at
     their target, so the targets are parented to the camera as well -- that
     keeps the whole rig in camera space. */
  scene.add(camera);
  const rig = new THREE.Group();
  camera.add(rig);
  const addLight = (color, intensity, pos) => {
    const l = new THREE.DirectionalLight(color, intensity);
    l.position.set(pos[0], pos[1], pos[2]);
    l.target.position.set(0, 0, -1);
    rig.add(l, l.target);
  };
  addLight(0xffffff, 1.55, [0.5, 0.7, 0.3]);    /* key, over the viewer's shoulder */
  addLight(0xcfdcf2, 0.75, [-0.85, -0.4, 0.25]); /* fill, opposite and cooler */
  addLight(0xffffff, 0.5, [0, 0.35, -2.5]);      /* rim, picks out the silhouette */
  /* soft omnidirectional base so nothing ever goes fully black */
  scene.add(new THREE.HemisphereLight(0xf4f6f9, 0x9aa1aa, 1.15));

  const root = new THREE.Group();
  scene.add(root);

  const mkMat = (o) => new THREE.MeshStandardMaterial(
    Object.assign({ roughness: 0.42, metalness: 0.03 }, o));

  const atoms = new THREE.InstancedMesh(
    new THREE.SphereGeometry(1, 20, 14), mkMat({}), MAXA);
  const bonds = new THREE.InstancedMesh(
    new THREE.CylinderGeometry(0.11, 0.11, 1, 10), mkMat({}), MAXB);
  const esp = new THREE.InstancedMesh(
    new THREE.SphereGeometry(1, 12, 9),
    mkMat({ transparent: true, opacity: 0.52, roughness: 0.6 }), n3);
  const pharm = new THREE.InstancedMesh(
    new THREE.SphereGeometry(1, 18, 12),
    mkMat({ transparent: true, opacity: 0.4, depthWrite: false }), MAXS);
  const shafts = new THREE.InstancedMesh(
    new THREE.CylinderGeometry(0.1, 0.1, 1, 8), mkMat({}), MAXR);
  const heads = new THREE.InstancedMesh(
    new THREE.ConeGeometry(0.2, 0.46, 10), mkMat({}), MAXR);

  const groups = { x1: new THREE.Group(), x3: new THREE.Group(), x4: new THREE.Group() };
  groups.x1.add(atoms, bonds);
  groups.x3.add(esp);
  groups.x4.add(pharm, shafts, heads);
  root.add(groups.x1, groups.x3, groups.x4);
  for (const m of [atoms, bonds, esp, pharm, shafts, heads]) m.frustumCulled = false;

  const controls = new TrackballControls(camera, renderer.domElement);
  controls.rotateSpeed = 3.2;
  controls.zoomSpeed = 1.1;
  controls.noPan = true;
  controls.staticMoving = false;
  controls.dynamicDampingFactor = 0.14;
  controls.minDistance = DIST * 0.3;
  controls.maxDistance = DIST * 3;
  /* TrackballControls binds A/S/D on the window by default; this page should
     not swallow those keys. */
  controls.keys = [];

  /* ---------- per-frame geometry ---------- */
  const M = new THREE.Matrix4(), Q = new THREE.Quaternion();
  const A = new THREE.Vector3(), B = new THREE.Vector3(), V = new THREE.Vector3();
  const UP = new THREE.Vector3(0, 1, 0);
  const C = new THREE.Color();
  const px = new Float32Array(MAXA * 3);

  function placeCylinder(mesh, i, x0, y0, z0, x1_, y1_, z1_, radiusScale) {
    A.set(x0, y0, z0); B.set(x1_, y1_, z1_);
    V.subVectors(B, A);
    const len = V.length() || 1e-6;
    Q.setFromUnitVectors(UP, V.divideScalar(len));
    M.compose(A.addVectors(A, V.multiplyScalar(len / 2)), Q,
              new THREE.Vector3(radiusScale, len, radiusScale));
    mesh.setMatrixAt(i, M);
  }

  function drawFrame(f) {
    /* --- atoms --- */
    const a0 = x1o[f], na = x1c[f];
    for (let i = 0; i < na; i++) {
      const e = x1e[a0 + i], st = elemStyle[e];
      const x = x1p[(a0 + i) * 3] / S, y = x1p[(a0 + i) * 3 + 1] / S, z = x1p[(a0 + i) * 3 + 2] / S;
      px[i * 3] = x; px[i * 3 + 1] = y; px[i * 3 + 2] = z;
      M.compose(A.set(x, y, z), Q.identity(), V.setScalar(st.r * 0.55));
      atoms.setMatrixAt(i, M);
      atoms.setColorAt(i, C.setHex(st.c));
    }
    atoms.count = na;
    atoms.instanceMatrix.needsUpdate = true;
    if (atoms.instanceColor) atoms.instanceColor.needsUpdate = true;

    /* --- bonds: distance test against covalent radii --- */
    let nb = 0;
    for (let i = 0; i < na && nb < MAXB - 1; i++) {
      const ri = elemCov[x1e[a0 + i]];
      for (let j = i + 1; j < na && nb < MAXB - 1; j++) {
        const dx = px[i * 3] - px[j * 3], dy = px[i * 3 + 1] - px[j * 3 + 1],
              dz = px[i * 3 + 2] - px[j * 3 + 2];
        const d2 = dx * dx + dy * dy + dz * dz;
        const cut = ri + elemCov[x1e[a0 + j]] + 0.4;
        if (d2 > cut * cut || d2 < 0.16) continue;
        const mx = (px[i * 3] + px[j * 3]) / 2, my = (px[i * 3 + 1] + px[j * 3 + 1]) / 2,
              mz = (px[i * 3 + 2] + px[j * 3 + 2]) / 2;
        placeCylinder(bonds, nb, px[i * 3], px[i * 3 + 1], px[i * 3 + 2], mx, my, mz, 1);
        bonds.setColorAt(nb, C.setHex(elemStyle[x1e[a0 + i]].c));
        nb++;
        placeCylinder(bonds, nb, mx, my, mz, px[j * 3], px[j * 3 + 1], px[j * 3 + 2], 1);
        bonds.setColorAt(nb, C.setHex(elemStyle[x1e[a0 + j]].c));
        nb++;
      }
    }
    bonds.count = nb;
    bonds.instanceMatrix.needsUpdate = true;
    if (bonds.instanceColor) bonds.instanceColor.needsUpdate = true;

    /* --- ESP surface --- */
    const e0 = f * n3;
    for (let i = 0; i < n3; i++) {
      M.compose(A.set(x3p[(e0 + i) * 3] / S, x3p[(e0 + i) * 3 + 1] / S, x3p[(e0 + i) * 3 + 2] / S),
                Q.identity(), V.setScalar(0.22));
      esp.setMatrixAt(i, M);
      esp.setColorAt(i, espColor(x3e[e0 + i] / 127, C));
    }
    esp.instanceMatrix.needsUpdate = true;
    if (esp.instanceColor) esp.instanceColor.needsUpdate = true;

    /* --- pharmacophore anchors --- */
    const p0 = s4o[f], ns = s4c[f];
    for (let i = 0; i < ns; i++) {
      M.compose(A.set(s4p[(p0 + i) * 3] / S, s4p[(p0 + i) * 3 + 1] / S, s4p[(p0 + i) * 3 + 2] / S),
                Q.identity(), V.setScalar(Math.max(0.3, s4rad[p0 + i] / 50)));
      pharm.setMatrixAt(i, M);
      pharm.setColorAt(i, palette[s4col[p0 + i]]);
    }
    pharm.count = ns;
    pharm.instanceMatrix.needsUpdate = true;
    if (pharm.instanceColor) pharm.instanceColor.needsUpdate = true;

    /* --- pharmacophore vectors --- */
    const r0 = a4o[f], nr = a4c[f];
    for (let i = 0; i < nr; i++) {
      const k = r0 + i;
      const ox = a4p[k * 3] / S, oy = a4p[k * 3 + 1] / S, oz = a4p[k * 3 + 2] / S;
      const L = a4len[k] / 20;
      const dx = a4d[k * 3] / 127, dy = a4d[k * 3 + 1] / 127, dz = a4d[k * 3 + 2] / 127;
      const sl = Math.max(0.01, L - 0.46);
      placeCylinder(shafts, i, ox, oy, oz, ox + dx * sl, oy + dy * sl, oz + dz * sl, 1);
      shafts.setColorAt(i, palette[a4col[k]]);
      V.set(dx, dy, dz).normalize();
      Q.setFromUnitVectors(UP, V);
      M.compose(A.set(ox + dx * (L - 0.23), oy + dy * (L - 0.23), oz + dz * (L - 0.23)),
                Q, B.set(1, 1, 1));
      heads.setMatrixAt(i, M);
      heads.setColorAt(i, palette[a4col[k]]);
    }
    shafts.count = nr; heads.count = nr;
    shafts.instanceMatrix.needsUpdate = true; heads.instanceMatrix.needsUpdate = true;
    if (shafts.instanceColor) shafts.instanceColor.needsUpdate = true;
    if (heads.instanceColor) heads.instanceColor.needsUpdate = true;
  }

  /* ---------- UI ---------- */
  const ui = document.createElement('div');
  ui.className = 'tv-ui';
  ui.innerHTML = `
    <div class="tv-row tv-transport">
      <button class="tv-btn tv-play" type="button" aria-label="Play">
        <span class="tv-ico-play" aria-hidden="true"></span><span class="tv-play-label">Play</span>
      </button>
      <input class="tv-slider" type="range" min="0" max="${nF - 1}" value="0" step="1"
             aria-label="Diffusion step">
      <span class="tv-step" aria-live="off">step 1</span>
    </div>
    <div class="tv-row tv-toggles">
      <label><input type="checkbox" data-layer="x1" checked> Molecule</label>
      <label><input type="checkbox" data-layer="x3" checked> ESP surface</label>
      <label><input type="checkbox" data-layer="x4" checked> Pharmacophores</label>
      <label><input type="checkbox" data-spin> Spin</label>
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
        speedSel = $('.tv-speed');

  let frame = 0, playing = false, spin = false, speed = 1, last = 0, acc = 0;
  const FPS = 28;

  /* The readout counts the denoising schedule, not the dumped frames: the
     first frame is step 1 and the last is DISPLAY_STEPS, with the dumped
     frames spread linearly in between. */
  const DISPLAY_STEPS = opts.displaySteps || 400;
  const srcLast = (data.srcFrames ? data.srcFrames - 1 : nF - 1) || 1;
  const stepLabel = (f) => {
    const t = data.frameIndex ? data.frameIndex[f] : f;
    return Math.round(1 + (t / srcLast) * (DISPLAY_STEPS - 1));
  };

  function setFrame(f, fromSlider) {
    frame = Math.max(0, Math.min(nF - 1, f | 0));
    drawFrame(frame);
    if (!fromSlider) slider.value = frame;
    stepOut.textContent = `step ${stepLabel(frame)} / ${DISPLAY_STEPS}`;
  }
  function setPlaying(p) {
    playing = p;
    playBtn.classList.toggle('is-playing', p);
    playBtn.setAttribute('aria-label', p ? 'Pause' : 'Play');
    $('.tv-play-label').textContent = p ? 'Pause' : 'Play';
  }

  playBtn.addEventListener('click', () => {
    if (!playing && frame >= nF - 1) setFrame(0);
    setPlaying(!playing);
  });
  slider.addEventListener('input', () => { setPlaying(false); setFrame(+slider.value, true); });
  speedSel.addEventListener('change', () => { speed = +speedSel.value; });
  ui.querySelectorAll('[data-layer]').forEach(cb => cb.addEventListener('change', () => {
    groups[cb.dataset.layer].visible = cb.checked;
  }));
  $('[data-spin]').addEventListener('change', (e) => { spin = e.target.checked; });
  $('.tv-reset').addEventListener('click', () => {
    controls.reset(); camera.position.set(0, 0, DIST); root.rotation.set(0, 0, 0);
  });
  box.addEventListener('keydown', (e) => {
    if (e.key === ' ') { e.preventDefault(); playBtn.click(); }
    else if (e.key === 'ArrowRight') { setPlaying(false); setFrame(frame + 1); }
    else if (e.key === 'ArrowLeft') { setPlaying(false); setFrame(frame - 1); }
  });
  box.tabIndex = 0;

  /* ---------- resize + loop ---------- */
  function resize() {
    const w = host.clientWidth, h = host.clientHeight;
    if (!w || !h) return;
    renderer.setSize(w, h, false);
    camera.aspect = w / h;
    camera.updateProjectionMatrix();
    controls.handleResize();   /* TrackballControls caches the element rect */
  }
  const ro = new ResizeObserver(resize); ro.observe(host);
  resize();

  let visible = true;
  const vo = new IntersectionObserver((es) => { visible = es[0].isIntersecting; },
    { threshold: 0.01 });
  vo.observe(host);

  let alive = true;
  function tick(now) {
    if (!alive) return;
    requestAnimationFrame(tick);
    if (!visible) { last = now; return; }
    const dt = Math.min(0.1, (now - last) / 1000) || 0; last = now;
    if (playing) {
      acc += dt * FPS * speed;
      if (acc >= 1) {
        const adv = Math.floor(acc); acc -= adv;
        if (frame + adv >= nF - 1) { setFrame(nF - 1); setPlaying(false); }
        else setFrame(frame + adv);
      }
    }
    if (spin) root.rotation.y += dt * 0.35;
    controls.update();
    renderer.render(scene, camera);
  }

  setFrame(0);
  requestAnimationFrame((t) => { last = t; tick(t); });

  const reduced = matchMedia('(prefers-reduced-motion: reduce)').matches;
  if (!reduced && opts.autoplay !== false) {
    const io = new IntersectionObserver((es) => {
      if (es[0].isIntersecting) { setPlaying(true); io.disconnect(); }
    }, { threshold: 0.4 });
    io.observe(host);
  }

  function dispose() {
    alive = false;
    ro.disconnect(); vo.disconnect();
    controls.dispose();
    for (const m of [atoms, bonds, esp, pharm, shafts, heads]) {
      m.geometry.dispose(); m.material.dispose(); m.dispose();
    }
    renderer.dispose();
    renderer.domElement.remove();
    ui.remove();
  }

  return { setFrame, setPlaying, dispose };
}
