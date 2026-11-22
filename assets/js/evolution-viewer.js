/* Ligand evolution viewer -----------------------------------------------------
   One scene, six stages.  Rather than three unrelated trajectories behind
   three tabs, this walks a single story: crystallographic fragments bound to
   a protein, two of them chosen and their interaction profiles extracted, a
   molecule generated to merge those profiles, a third fragment brought in, a
   second generation,
   and the final molecule over all three.

   Every structure lives in the 8CNX crystal frame, so nothing is re-centred
   between stages -- moving from stage to stage only changes what is visible
   and how opaque it is, never where anything sits.  That is the whole point:
   the overlays line up because they really do line up.

   Data: tools/build_evolution.py -> evo-base.js, evo-traj-1.js, evo-traj-2.js
---------------------------------------------------------------------------- */
import * as THREE from 'three';
import { TrackballControls } from 'three/addons/controls/TrackballControls.js';
import { createActor } from './mol-core.js';
import { buildRibbon } from './ribbon.js';
import { principalAxes, leastCrowded, uprightFor } from './viewpoint.js';

const FOV = 38;
const DISPLAY_STEPS = 400;
const FPS = 28;

/* Carbon colours, so molecules sharing a view can be told apart the way a
   medicinal-chemistry figure does it.  Heteroatoms keep their CPK colours in
   every structure, which is what makes the overlays comparable. */
const CARBON = {
  frag8:   0x2f9e8f,
  frag10:  0xd4902f,
  frag12:  0x9478cf,
  merged1: 0x4c5057,
  merged2: 0x4c5057,
  t1:      0x4c5057,
  t2:      0x4c5057,
};

/* lead = what the stage is about, drawn solid and allowed to show its profile.
   refs = context, drawn as ghosts with no profile of their own. */
const STAGES = [
  { label: 'Bound fragments', lead: ['frag8', 'frag10', 'frag12'], refs: [],
    profiles: false, wide: true,
    note: 'A fragment screen against the D68 3C protease: three fragments soak into the same pocket, each making only a few contacts.' },
  { label: 'Pick 1 and 2', lead: ['frag8', 'frag10'], refs: [],
    profiles: true,
    note: 'Fragments 1 and 2 occupy neighbouring subsites. Fragment 3 is set aside. The interaction profiles of the two — ESP surface and pharmacophores — are extracted and combined into one target.' },
  { label: 'Merge round 1', traj: 't1', lead: ['t1'], refs: ['frag8', 'frag10'],
    profiles: false,
    note: 'ShEPhERD-2 denoises the molecule and its interaction profile together, conditioned on the combined profile of fragments 1 and 2. Only the atoms are drawn here — the profile is still noise for most of the run.' },
  { label: 'Pick 3', lead: ['merged1'], refs: ['frag12'],
    profiles: false,
    note: 'Fragment 3 binds an adjacent subsite the round-one molecule does not reach. Its profile is added to the target.' },
  { label: 'Merge round 2', traj: 't2', lead: ['t2'], refs: ['merged1', 'frag12'],
    profiles: false,
    note: 'A second generation, conditioned on the round-one molecule plus fragment 3 — the merged product goes back into the model as input.' },
  { label: 'Final molecule', lead: ['merged2'], refs: ['frag8', 'frag10', 'frag12'],
    profiles: false,
    note: 'The round-two product over all three starting fragments: one molecule covering three fragment subsites.' },
];

const REF_OPACITY = 0.38;
const EASE = 6.5;              /* opacity/camera convergence rate, per second */

export function mount(box, base, opts = {}) {
  const reduced = matchMedia('(prefers-reduced-motion: reduce)').matches;
  const ZERO = new THREE.Vector3(0, 0, 0);
  const stage = box.querySelector('.tv-stage') || box.appendChild(
    Object.assign(document.createElement('div'), { className: 'tv-stage' }));
  stage.querySelectorAll('.tv-fallback').forEach(n => n.remove());
  const host = stage;

  /* ---------- scene ---------- */
  const scene = new THREE.Scene();
  const camera = new THREE.PerspectiveCamera(FOV, 1, 0.5, 4000);
  const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
  renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
  renderer.outputColorSpace = THREE.SRGBColorSpace;
  host.appendChild(renderer.domElement);
  renderer.domElement.classList.add('tv-canvas');

  /* Lighting rides on the camera, so whichever way the scene is turned the
     face you are looking at is the face that is lit.  Directional lights aim
     at a target, so the targets are parented to the camera too. */
  scene.add(camera);
  const rig = new THREE.Group();
  camera.add(rig);
  const addLight = (color, intensity, pos) => {
    const l = new THREE.DirectionalLight(color, intensity);
    l.position.set(pos[0], pos[1], pos[2]);
    l.target.position.set(0, 0, -1);
    rig.add(l, l.target);
  };
  addLight(0xffffff, 1.55, [0.5, 0.7, 0.3]);
  addLight(0xcfdcf2, 0.75, [-0.85, -0.4, 0.25]);
  addLight(0xffffff, 0.5, [0, 0.35, -2.5]);
  scene.add(new THREE.HemisphereLight(0xf4f6f9, 0x9aa1aa, 1.15));

  const root = new THREE.Group();
  scene.add(root);

  /* ---------- protein ---------- */
  /* farA is generous because only the outer wall of the tube is drawn now:
     the far side of the fold has to carry itself on one layer, not two */
  const ribbonGeom = buildRibbon(base.protein, { nearR: 12, farR: 34, farA: 0.58 });
  const ribbonMat = new THREE.MeshStandardMaterial({
    color: 0xcfd4dc, roughness: 0.62, metalness: 0.0,
    transparent: true, opacity: 0.3, depthWrite: false,
    /* the tube is closed and wound outward, so only its outer wall is drawn:
       a DoubleSide translucent tube blends its own back wall over its front
       and every crossing of the chain reads twice as dark as it should */
    side: THREE.FrontSide, vertexColors: true,
  });
  const ribbon = new THREE.Mesh(ribbonGeom, ribbonMat);
  ribbon.renderOrder = 10;           /* after the molecules, so it tints them */
  ribbon.frustumCulled = false;
  root.add(ribbon);

  /* where the protein sits as a whole, for the wide opening shot */
  const proteinCenter = ribbonGeom.boundingSphere.center.clone();
  const proteinRadius = ribbonGeom.boundingSphere.radius;

  /* ---------- actors ---------- */
  const shared = { elements: base.elements, colors: base.colors };
  const actors = {};
  for (const [key, data] of Object.entries(base.structures)) {
    const a = createActor(Object.assign({ posScale: base.posScale }, data),
                          Object.assign({ carbon: CARBON[key] }, shared));
    a.group.visible = false;
    root.add(a.group);
    actors[key] = a;
  }

  /* A trajectory is several hundred kilobytes, so it is fetched the first time
     its stage is opened rather than on page load. */
  const trajPending = {};
  function ensureTraj(id) {
    if (actors[id]) return Promise.resolve(actors[id]);
    if (trajPending[id]) return trajPending[id];
    trajPending[id] = Promise.resolve(opts.loadTraj(id)).then((data) => {
      const a = createActor(Object.assign({ posScale: base.posScale }, data),
                            Object.assign({ carbon: CARBON[id] }, shared));
      a.group.visible = false;
      a.setShowH(showH);
      root.add(a.group);
      actors[id] = a;
      return a;
    });
    return trajPending[id];
  }

  /* ---------- camera framing ---------- */
  const fitDist = (r) => (r * 1.15) / Math.tan((FOV / 2) * Math.PI / 180);
  let pocketR = 1;
  for (const a of Object.values(actors)) pocketR = Math.max(pocketR, a.radius(0, true));
  const POCKET_DIST = fitDist(pocketR);
  /* The opening shot frames the whole fold with the pocket sitting low and
     centre, which is where the site lands when the camera looks at the
     protein's own centre from outside.  Cropping below the bounding radius
     lets the fold fill the frame rather than float in it. */
  const WIDE_DIST = fitDist(proteinRadius * 0.88);
  /* not quite the protein's centre: pulled a little toward the pocket, which
     lifts the fragments off the bottom edge without losing the whole fold */
  const WIDE_TARGET = proteinCenter.clone().multiplyScalar(0.72);

  /* Opening direction: look into the pocket from outside the protein, then
     nudge within a cap of that to the angle that piles the final molecule's
     heavy atoms up the least.  Starting from the outward normal is what keeps
     the camera from opening buried inside the fold. */
  const HOME_DIR = proteinCenter.clone().negate().normalize();
  const HOME_UP = new THREE.Vector3(0, 1, 0);
  {
    const pts = [];
    const d = actors.merged2.data, S = base.posScale;
    const el = base.elements;
    /* heavy atoms of the final molecule, straight out of the payload */
    const elemIdx = atobBytes(d.x1.elem);
    const posArr = atobInt16(d.x1.pos);
    for (let i = 0; i < elemIdx.length; i++) {
      if (el[elemIdx[i]] === 'H') continue;
      pts.push([posArr[i * 3] / S, posArr[i * 3 + 1] / S, posArr[i * 3 + 2] / S]);
    }
    const ax = principalAxes(pts);
    if (ax) {
      const dir = leastCrowded(ax.pts, HOME_DIR, THREE.MathUtils.degToRad(38), ax);
      HOME_DIR.copy(dir);
      HOME_UP.copy(uprightFor(ax.pts, dir));
    }
  }

  const target = new THREE.Vector3(0, 0, 0);        /* where the camera looks */
  const wantTarget = new THREE.Vector3(0, 0, 0);
  let wantDist = POCKET_DIST;
  camera.up.copy(HOME_UP);
  camera.position.copy(HOME_DIR).multiplyScalar(POCKET_DIST);
  camera.lookAt(target);

  const controls = new TrackballControls(camera, renderer.domElement);
  controls.rotateSpeed = 3.2;
  controls.zoomSpeed = 1.1;
  controls.noPan = true;
  controls.staticMoving = false;
  controls.dynamicDampingFactor = 0.14;
  controls.minDistance = POCKET_DIST * 0.25;
  controls.maxDistance = WIDE_DIST * 1.8;
  controls.keys = [];                /* do not swallow A/S/D on the page */
  controls.target = target;

  /* ---------- stage state ---------- */
  let current = -1;
  let showH = false;
  let showProtein = true;
  let showProfiles = false;
  let spin = false;
  let speed = 1;
  let camAnim = false;
  /* each animated stage keeps its own position in its own trajectory */
  const play = { t1: { frame: 0, playing: false, acc: 0 },
                 t2: { frame: 0, playing: false, acc: 0 } };

  /* Target opacities, eased every tick so stages dissolve into each other
     instead of cutting. */
  const want = {};       /* key -> {mol, prof} */
  const have = {};
  for (const k of ['frag8', 'frag10', 'frag12', 'merged1', 'merged2', 't1', 't2']) {
    want[k] = { mol: 0, prof: 0 };
    have[k] = { mol: 0, prof: 0 };
  }
  let wantRibbon = 0.3, haveRibbon = 0.3;

  /* Only the fragments and the finished molecules carry an interaction
     profile.  The generation trajectories are atoms alone. */
  const hasProfile = (k) => !!(actors[k] && actors[k].hasProfile);
  const stageHasProfile = (st) => !!st && st.lead.some(hasProfile);

  function refreshTargets() {
    const st = STAGES[current];
    for (const k of Object.keys(want)) want[k].mol = want[k].prof = 0;
    for (const k of st.lead) {
      want[k].mol = 1;
      /* a trajectory carries no profile of its own, so there is nothing to
         raise here however the toggle is set */
      want[k].prof = showProfiles && hasProfile(k) ? 1 : 0;
    }
    for (const k of st.refs) want[k].mol = REF_OPACITY;
    /* The wide shot is the one place the whole fold should read; up close the
       ribbon is backdrop and must not wash the molecules out. */
    wantRibbon = showProtein ? (st.wide ? 0.55 : 0.38) : 0;
  }

  function applyOpacity(k, o) {
    const a = actors[k];
    if (!a) return;
    a.group.visible = o.mol > 0.004 || o.prof > 0.004;
    a.setLayerOpacity('x1', o.mol);
    a.setLayerOpacity('x3', o.prof);
    a.setLayerOpacity('x4', o.prof);
  }

  function goTo(i, { snap = false } = {}) {
    if (i === current) return;
    current = i;
    const st = STAGES[i];
    showProfiles = st.profiles;
    refreshTargets();
    wantDist = st.wide ? WIDE_DIST : POCKET_DIST;
    wantTarget.copy(st.wide ? WIDE_TARGET : ZERO);
    /* Only a stage change may move the camera.  Once it has arrived the tween
       switches itself off, so a viewer's own zoom is not dragged back. */
    camAnim = true;
    if (st.traj) {
      ensureTraj(st.traj).then((a) => {
        a.setFrame(play[st.traj].frame);
        refreshTargets();
        syncUI();
        /* opening an animation stage for the first time starts it running */
        if (!reduced && !play[st.traj].started) {
          play[st.traj].started = true;
          play[st.traj].playing = true;
          syncUI();
        }
      }).catch((e) => { console.error(e); syncUI(); });
    }
    for (const p of Object.values(play)) p.playing = false;
    if (snap) {
      for (const k of Object.keys(want)) { have[k] = Object.assign({}, want[k]); applyOpacity(k, have[k]); }
      haveRibbon = wantRibbon;
      target.copy(wantTarget);
      camera.position.copy(HOME_DIR).multiplyScalar(wantDist).add(target);
      camera.up.copy(HOME_UP);
      camAnim = false;
    }
    syncUI();
  }

  /* ---------- UI ---------- */
  const ui = document.createElement('div');
  ui.className = 'evo-ui';
  ui.innerHTML = `
    <ol class="evo-timeline">${STAGES.map((s, i) => `
      <li class="evo-stage${s.traj ? ' is-anim' : ''}" data-stage="${i}">
        <button class="evo-chip" type="button">
          <span class="evo-num">${i + 1}</span>
          <span class="evo-name">${s.label}</span>
        </button>
        ${s.traj ? `
        <div class="evo-track">
          <button class="evo-play" type="button" aria-label="Play">
            <span class="evo-ico" aria-hidden="true"></span>
          </button>
          <input class="evo-slider" type="range" min="0" max="100" value="0" step="1"
                 aria-label="Denoising step">
          <span class="evo-step">1 / ${DISPLAY_STEPS}</span>
          <span class="evo-mini" aria-hidden="true"><i></i></span>
        </div>` : ''}
      </li>`).join('')}</ol>
    <p class="evo-note" id="evo-note"></p>
    <div class="evo-controls">
      <label><input type="checkbox" data-opt="protein" checked> Protein</label>
      <label><input type="checkbox" data-opt="profiles"> Interaction profile</label>
      <label><input type="checkbox" data-opt="hydrogens"> Hydrogens</label>
      <label><input type="checkbox" data-opt="spin"> Spin</label>
      <select class="evo-speed" aria-label="Playback speed">
        <option value="0.5">0.5&times;</option>
        <option value="1" selected>1&times;</option>
        <option value="2">2&times;</option>
      </select>
      <button class="evo-btn" type="button" data-act="reset">Reset view</button>
    </div>`;
  box.appendChild(ui);

  const noteEl = ui.querySelector('.evo-note');
  const stageEls = [...ui.querySelectorAll('.evo-stage')];
  const optEl = (n) => ui.querySelector(`[data-opt="${n}"]`);

  stageEls.forEach((el, i) => {
    el.querySelector('.evo-chip').addEventListener('click', () => goTo(i));
    const st = STAGES[i];
    if (!st.traj) return;
    const sl = el.querySelector('.evo-slider');
    const pb = el.querySelector('.evo-play');
    sl.addEventListener('input', () => {
      if (current !== i) goTo(i);
      const a = actors[st.traj];
      if (!a) return;
      play[st.traj].playing = false;
      play[st.traj].frame = Math.round(+sl.value / 100 * (a.nframes - 1));
      a.setFrame(play[st.traj].frame);
      syncUI();
    });
    pb.addEventListener('click', () => {
      if (current !== i) { goTo(i); return; }
      const a = actors[st.traj], p = play[st.traj];
      if (!a) return;
      if (!p.playing && p.frame >= a.nframes - 1) { p.frame = 0; a.setFrame(0); }
      p.playing = !p.playing;
      p.started = true;
      syncUI();
    });
  });

  optEl('protein').addEventListener('change', (e) => {
    showProtein = e.target.checked; refreshTargets();
  });
  optEl('profiles').addEventListener('change', (e) => {
    showProfiles = e.target.checked; refreshTargets();
  });
  optEl('hydrogens').addEventListener('change', (e) => {
    showH = e.target.checked;
    for (const a of Object.values(actors)) a.setShowH(showH);
  });
  optEl('spin').addEventListener('change', (e) => { spin = e.target.checked; });
  ui.querySelector('.evo-speed').addEventListener('change', (e) => { speed = +e.target.value; });
  ui.querySelector('[data-act="reset"]').addEventListener('click', () => {
    controls.reset();
    target.copy(wantTarget);
    camera.position.copy(HOME_DIR).multiplyScalar(wantDist).add(target);
    camera.up.copy(HOME_UP);
    camera.lookAt(target);
  });

  box.tabIndex = 0;
  box.addEventListener('keydown', (e) => {
    const st = STAGES[current];
    if (e.key === 'ArrowRight' && !st.traj) { goTo(Math.min(STAGES.length - 1, current + 1)); }
    else if (e.key === 'ArrowLeft' && !st.traj) { goTo(Math.max(0, current - 1)); }
    else if (e.key === ' ' && st.traj) {
      e.preventDefault();
      stageEls[current].querySelector('.evo-play').click();
    } else if (st.traj && (e.key === 'ArrowRight' || e.key === 'ArrowLeft')) {
      const a = actors[st.traj], p = play[st.traj];
      if (!a) return;
      p.playing = false;
      p.frame = Math.max(0, Math.min(a.nframes - 1, p.frame + (e.key === 'ArrowRight' ? 1 : -1)));
      a.setFrame(p.frame);
      syncUI();
    } else return;
    e.preventDefault();
  });

  /* The readout counts the denoising schedule, not the frames that were kept:
     the first is step 1 and the last is DISPLAY_STEPS, spread linearly. */
  function stepLabel(a, f) {
    const srcLast = (a.data.srcFrames ? a.data.srcFrames - 1 : a.nframes - 1) || 1;
    const t = a.data.frameIndex ? a.data.frameIndex[f] : f;
    return Math.round(1 + (t / srcLast) * (DISPLAY_STEPS - 1));
  }

  function syncUI() {
    const st = STAGES[current];
    noteEl.textContent = st ? st.note : '';
    stageEls.forEach((el, i) => {
      el.classList.toggle('is-active', i === current);
      el.setAttribute('aria-current', i === current ? 'step' : 'false');
      const s = STAGES[i];
      if (!s.traj) return;
      const a = actors[s.traj], p = play[s.traj];
      const frac = a ? p.frame / Math.max(1, a.nframes - 1) : 0;
      el.querySelector('.evo-slider').value = Math.round(frac * 100);
      el.querySelector('.evo-mini').style.setProperty('--p', `${frac * 100}%`);
      el.querySelector('.evo-step').textContent =
        a ? `${stepLabel(a, p.frame)} / ${DISPLAY_STEPS}` : 'loading…';
      el.classList.toggle('is-playing', !!p.playing);
      el.querySelector('.evo-play')
        .setAttribute('aria-label', p.playing ? 'Pause' : 'Play');
    });
    /* nothing to toggle on a stage whose molecule has no profile: the control
       greys out rather than sitting there doing nothing when clicked */
    const canProf = stageHasProfile(st);
    const pe = optEl('profiles');
    pe.checked = showProfiles && canProf;
    pe.disabled = !canProf;
    pe.closest('label').classList.toggle('is-off', !canProf);
  }

  /* ---------- resize + loop ---------- */
  function resize() {
    const w = host.clientWidth, h = host.clientHeight;
    if (!w || !h) return;
    renderer.setSize(w, h, false);
    camera.aspect = w / h;
    camera.updateProjectionMatrix();
    controls.handleResize();          /* TrackballControls caches the rect */
  }
  const ro = new ResizeObserver(resize); ro.observe(host);
  resize();

  let onScreen = true;
  const vo = new IntersectionObserver((es) => { onScreen = es[0].isIntersecting; },
    { threshold: 0.01 });
  vo.observe(host);

  const tmp = new THREE.Vector3();
  let alive = true, last = 0;

  function tick(now) {
    if (!alive) return;
    requestAnimationFrame(tick);
    if (!onScreen) { last = now; return; }
    const dt = Math.min(0.1, (now - last) / 1000) || 0; last = now;
    const k = 1 - Math.exp(-EASE * dt);     /* frame-rate independent easing */

    /* playback */
    const st = STAGES[current];
    if (st && st.traj && actors[st.traj]) {
      const a = actors[st.traj], p = play[st.traj];
      if (p.playing) {
        p.acc += dt * FPS * speed;
        if (p.acc >= 1) {
          const adv = Math.floor(p.acc); p.acc -= adv;
          if (p.frame + adv >= a.nframes - 1) { p.frame = a.nframes - 1; p.playing = false; }
          else p.frame += adv;
          a.setFrame(p.frame);
          syncUI();
        }
      }
    }

    /* opacity */
    for (const key of Object.keys(want)) {
      const w = want[key], h = have[key];
      if (Math.abs(w.mol - h.mol) < 0.002 && Math.abs(w.prof - h.prof) < 0.002) {
        if (h.mol !== w.mol || h.prof !== w.prof) {
          h.mol = w.mol; h.prof = w.prof; applyOpacity(key, h);
        }
        continue;
      }
      h.mol += (w.mol - h.mol) * k;
      h.prof += (w.prof - h.prof) * k;
      applyOpacity(key, h);
    }
    if (Math.abs(wantRibbon - haveRibbon) > 0.002) {
      haveRibbon += (wantRibbon - haveRibbon) * k;
      ribbonMat.opacity = haveRibbon;
      ribbon.visible = haveRibbon > 0.004;
    }

    /* camera: ease the orbit centre and the distance, never the orientation --
       whatever the viewer has rotated to is theirs to keep across stages */
    if (camAnim) {
      tmp.subVectors(camera.position, target);
      const dist = tmp.length() || 1;
      const nd = dist + (wantDist - dist) * k * 0.8;
      target.lerp(wantTarget, k * 0.8);
      camera.position.copy(target).addScaledVector(tmp.divideScalar(dist), nd);
      if (target.distanceToSquared(wantTarget) < 1e-4 && Math.abs(nd - wantDist) < 1e-2)
        camAnim = false;
    }

    if (spin) {
      tmp.subVectors(camera.position, target);
      tmp.applyAxisAngle(camera.up, dt * 0.3);
      camera.position.copy(target).add(tmp);
    }

    controls.update();
    renderer.render(scene, camera);
  }

  goTo(0, { snap: true });
  requestAnimationFrame((t) => { last = t; tick(t); });

  function dispose() {
    alive = false;
    ro.disconnect(); vo.disconnect();
    controls.dispose();
    for (const a of Object.values(actors)) a.dispose();
    ribbonGeom.dispose(); ribbonMat.dispose();
    renderer.dispose();
    renderer.domElement.remove();
    ui.remove();
  }

  return { goTo, dispose };
}

/* small base64 helpers, used once at startup to pick the opening viewpoint */
function atobBytes(s) {
  const bin = atob(s), out = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) out[i] = bin.charCodeAt(i);
  return out;
}
function atobInt16(s) {
  const b = atobBytes(s);
  return new Int16Array(b.buffer, 0, b.length / 2);
}
