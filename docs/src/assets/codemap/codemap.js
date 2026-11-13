// ThetaL's code map (configured in the block below; data from docs/codemap.jl): the schemes on the left,
// then the options marked with the schemes they combine
// with, the directories of src/
// on the right (folds down to files, one dot per top-level function), and a curve from each setting to the visible
// rows whose functions it runs beyond its scheme (a scheme: beyond the base case, which is joined to everything it
// runs), its width growing with their number. Hovering or clicking a setting
// colours the functions it runs; hovering a dot names the function and the settings that run it. The list below
// the graph spells out the selection.
(function () {
  // ---- configuration for ThetaL ----
  const SCHEMEGROUP = "Scheme";                                    // the group whose settings the options combine with
  const GROUPS = [SCHEMEGROUP, "Velocity", "Integrator", "Driver", "Curves"];
  // modules, each with its directories; all start closed
  const MODULES = [
    { id: "ThetaL", about: "boundary evolution in θ–L variables",
      dirs: ["src", "curves", "stepping", "problem", "integrator", "drivers"] },
  ];
  const DIRNAME = { src: "src (periodic grid functions)" };
  const SCHEMEMARK = { etdrk2: "2", expeuler: "1" };
  const OPTIONS_LABEL = "Options (2 1: the schemes they combine with)";
  // ---- end of configuration ----
  const dirname = dir => DIRNAME[dir] || dir.split("/").pop();
  const SVG = "http://www.w3.org/2000/svg";

  function el(tag, cls, text) {
    const e = document.createElement(tag);
    if (cls) e.className = cls;
    if (text !== undefined) e.textContent = text;
    return e;
  }

  function setup(root, d) {
    const S = d.settings, F = d.functions, NS = S.length;
    const base = S.find(s => s.base);
    const sidx = Object.fromEntries(S.map((s, i) => [s.id, i]));
    const gcol = s => `var(--cm-g${GROUPS.indexOf(s.group)})`;
    const runs = S.map(() => new Set());                     // setting → indices of the functions it runs
    F.forEach((f, i) => f.used.forEach(id => runs[sidx[id]].add(i)));
    // each setting is compared with what it changes: an option with its scheme (`parent`), a scheme with the
    // base case; the base case shows everything it runs
    const ref = s => s.base ? null : S[sidx[s.parent || base.id]];
    const extra = s => s.base ? runs[sidx[s.id]] :
      new Set([...runs[sidx[s.id]]].filter(i => !runs[sidx[ref(s).id]].has(i)));
    const extras = S.map(extra);
    const beyond = s => s.base ? "" : ` beyond the ${ref(s).label.toLowerCase()}`;

    root.textContent = "";
    const graph = el("div", "cm-graph");
    const svg = document.createElementNS(SVG, "svg");
    svg.setAttribute("class", "cm-edges");
    const left = el("div", "cm-left"), right = el("div", "cm-right");
    graph.append(svg, left, right);
    const legend = el("div", "cm-legend");
    const detail = el("div", "cm-detail");
    const tip = el("div", "cm-tip");
    root.append(graph, legend, detail, tip);
    if (d.generated) root.append(el("p", "cm-hint cm-stamp", "Coverage of " + d.generated + "."));

    // settings: the schemes, then the options, which combine with the schemes marked on them; each option was run
    // with one scheme (`parent`) and is compared with it
    const buttons = [];
    const SCHEMES = S.filter(s => !s.parent);
    const schemeMark = SCHEMEMARK;
    function button(s, box) {
      const b = el("button", "cm-setting", s.label);
      b.type = "button";
      b.style.setProperty("--cm-c", gcol(s));
      b.dataset.id = s.id;
      b.title = s.about;
      if (s.parent) {
        const marks = el("span", "cm-schemes");
        SCHEMES.forEach(l => {
          const ok = (s.schemes || []).includes(l.id);
          const m = el("span", "cm-scheme" + (ok ? " ok" : "") + (l.id === s.parent ? " run" : ""), schemeMark[l.id] || l.label[0]);
          m.title = l.label + (l.id === s.parent ? ": run with it here" : ok ? ": combines" : ": does not combine");
          marks.append(m);
        });
        b.append(marks);
        b.title += `\nRun with the ${S[sidx[s.parent]].label.toLowerCase()}; combines with: ` +
          SCHEMES.filter(l => (s.schemes || []).includes(l.id)).map(l => l.label.toLowerCase()).join(", ") + ".";
      }
      buttons[sidx[s.id]] = b;
      box.append(b);
    }
    GROUPS.forEach(g => {
      const ss = S.filter(s => s.group === g);
      if (!ss.length) return;
      const box = el("div", "cm-group" + (g === SCHEMEGROUP ? " cm-schemes-group" : ""));
      box.append(el("div", "cm-group-label", g));
      ss.forEach(s => button(s, box));
      left.append(box);
      if (g === SCHEMEGROUP) left.append(el("div", "cm-group-label cm-options-label", OPTIONS_LABEL));
    });

    // files, by module and directory, one dot per function. Modules and their directories are folds: a closed
    // fold is one row (the curves end there), an open one shows what it holds.
    const fileRow = {}, fileFolds = {}, dots = [], folds = [];
    function fold(parent, name, about, fns, sub) {
      const head = el("button", "cm-mod" + (sub ? " cm-sub" : ""));
      head.type = "button";
      const count = el("span", "cm-mod-count");
      head.append(el("span", "cm-mod-name", name), el("span", "cm-mod-about", about ? " " + about : ""), count);
      const body = el("div", "cm-mod-body");
      const f = { head, body, count, fns, open: false };
      head.addEventListener("click", () => { f.open = !f.open; toggle(f); buildEdges(); show(); });
      parent.append(head, body);
      folds.push(f);
      toggle(f);
      return f;
    }
    function toggle(f) {
      f.head.classList.toggle("open", f.open);
      f.head.setAttribute("aria-expanded", f.open);
      f.body.hidden = !f.open;
    }
    const fnsOf = dirs => F.map((f, i) => i).filter(i => dirs.includes(d.files.find(x => x.path === F[i].file)?.dir));
    MODULES.forEach(m => {
      const mf = fold(right, m.id, m.about, fnsOf(m.dirs), false);
      m.dirs.forEach(dir => {
        const files = d.files.filter(f => f.dir === dir);
        if (!files.length) return;
        // a module of one directory shows its files directly
        const df = m.dirs.length > 1 ? fold(mf.body, dirname(dir), "", fnsOf([dir]), true) : null;
        files.forEach(file => {
          const row = el("div", "cm-file");
          const name = el("div", "cm-file-name", file.path.split("/").pop());
          name.title = "src/" + file.path;
          const box = el("div", "cm-dots");
          F.forEach((f, i) => {
            if (f.file !== file.path) return;
            const dot = el("span", "cm-dot");
            dot.dataset.i = i;
            dots[i] = dot;
            box.append(dot);
          });
          row.append(name, box);
          fileRow[file.path] = row;
          fileFolds[file.path] = df ? [mf, df] : [mf];
          (df || mf).body.append(row);
        });
      });
    });
    // the row a function shows on: its file, or the outermost closed fold around it
    const target = file => (fileFolds[file].find(f => !f.open) || {}).head || fileRow[file];

    // the curves: setting → visible row, width by the functions run there (beyond its scheme)
    let edges = [];                                           // {k, target, n, p}
    function buildEdges() {
      svg.textContent = "";
      edges = [];
      S.forEach((s, k) => {
        const per = new Map();
        extras[k].forEach(i => {
          const t = target(F[i].file);
          per.set(t, (per.get(t) || 0) + 1);
        });
        per.forEach((n, target) => {
          const p = document.createElementNS(SVG, "path");
          p.setAttribute("stroke-width", Math.min(1 + 1.2 * Math.sqrt(n), 7).toFixed(2));
          p.style.stroke = gcol(s);
          svg.append(p);
          edges.push({ k, target, n, p });
        });
      });
      layout();
    }
    function layout() {
      const g = graph.getBoundingClientRect();
      svg.setAttribute("width", g.width); svg.setAttribute("height", g.height);
      edges.forEach(e => {
        const b = buttons[e.k].getBoundingClientRect(), r = e.target.getBoundingClientRect();
        const a = (e.target.querySelector(".cm-mod-name, .cm-file-name") || e.target).getBoundingClientRect();
        const x0 = b.right - g.left, y0 = b.top + b.height / 2 - g.top;
        const x1 = r.left - g.left - 4, y1 = (a.top + a.bottom) / 2 - g.top;
        const c = (x1 - x0) / 2;
        e.p.setAttribute("d", `M${x0},${y0} C${x0 + c},${y0} ${x1 - c},${y1} ${x1},${y1}`);
      });
    }
    new ResizeObserver(layout).observe(graph);

    // state: a selected setting (click) and a previewed one (hover)
    let selected = null, preview = null;
    function show() {
      const k = preview ?? selected;
      root.classList.toggle("cm-active", k !== null);
      buttons.forEach((b, j) => b.classList.toggle("on", j === k));
      const scheme = k !== null && !S[k].parent ? S[k].id : null;
      buttons.forEach((b, j) => b.classList.toggle("na", scheme !== null && !!S[j].parent && !(S[j].schemes || []).includes(scheme)));
      edges.forEach(e => e.p.setAttribute("class", k === null ? "" : e.k === k ? "on" : "off"));
      F.forEach((f, i) => {
        const dot = dots[i];
        if (!dot) return;
        dot.style.background = ""; dot.style.borderColor = "";
        if (k === null) {
          const n = f.used.length;
          dot.className = "cm-dot " + (n === 0 ? "none" : n === NS ? "all" : "some");
          if (n > 0 && n < NS) dot.style.opacity = (0.3 + 0.7 * n / NS).toFixed(2);
          else dot.style.opacity = "";
        } else {
          dot.style.opacity = "";
          if (extras[k].has(i)) {
            dot.className = "cm-dot hit";
            dot.style.background = gcol(S[k]); dot.style.borderColor = gcol(S[k]);
          } else dot.className = "cm-dot " + (runs[k].has(i) ? "shared" : "off");
        }
      });
      folds.forEach(f => {
        const all = f.fns.length;
        const n = k === null ? 0 : f.fns.filter(i => extras[k].has(i)).length;
        f.count.textContent = k === null ? ` · ${all} functions` :
          ` · ${n} of ${all}` + (f.head.classList.contains("cm-sub") ? "" : beyond(S[k]));
        f.head.style.setProperty("--cm-c", k !== null && n ? gcol(S[k]) : "");
        f.head.classList.toggle("hit", k !== null && n > 0);
      });
      showDetail(k);
      showLegend(k);
    }

    function chips(list, box) {
      const byFile = {};
      list.forEach(i => (byFile[F[i].file] = byFile[F[i].file] || []).push(i));
      Object.keys(byFile).sort((a, b) => d.files.findIndex(f => f.path === a) - d.files.findIndex(f => f.path === b))
        .forEach(file => {
          const line = el("div", "cm-list-row");
          line.append(el("span", "cm-list-file", file));
          const names = el("span", "cm-list-names");
          byFile[file].sort((a, b) => F[a].line - F[b].line).forEach(i => {
            const c = el("code", "cm-chip", F[i].name);
            c.title = `src/${F[i].file}:${F[i].line}` + (F[i].doc ? "\n" + F[i].doc : "");
            names.append(c);
          });
          line.append(names);
          box.append(line);
        });
    }

    function showDetail(k) {
      detail.textContent = "";
      if (k === null) {
        const all = F.map((f, i) => i).filter(i => F[i].used.length === NS);
        const none = F.map((f, i) => i).filter(i => F[i].used.length === 0);
        detail.append(el("p", "cm-hint", "Hover over a setting to preview it; click one to keep it selected."));
        const a = el("details");
        a.append(el("summary", "", `Run by every setting (${all.length})`));
        chips(all, a);
        const b = el("details");
        b.append(el("summary", "", `Run by none of these settings (${none.length}): other options, display, input and output`));
        chips(none, b);
        detail.append(a, b);
        return;
      }
      const s = S[k];
      const head = el("p");
      const t = el("strong", "", s.label);
      t.style.color = gcol(s);
      head.append(t, document.createTextNode(` (${s.about}). Runs ${runs[k].size} of ${F.length} functions` +
        (s.base ? "." : `, ${extras[k].size} of them${beyond(s)}.`)));
      detail.append(head);
      if (!extras[k].size)
        detail.append(el("p", "cm-hint", `Nothing${beyond(s)}: the same functions run, with other data.`));
      chips([...extras[k]], detail);
    }

    function showLegend(k) {
      legend.textContent = "";
      const item = (cls, text, colour) => {
        const it = el("span", "cm-legend-item");
        const dot = el("span", "cm-dot " + cls);
        if (colour) { dot.style.background = colour; dot.style.borderColor = colour; }
        it.append(dot, document.createTextNode(text));
        legend.append(it);
      };
      if (k === null) {
        item("all", "run by every setting");
        item("some", "run by some (stronger: more of them)");
        item("none", "run by none");
      } else if (S[k].base) {
        item("hit", "run by the base case", gcol(S[k]));
        item("off", "not run");
      } else {
        item("hit", `run by ${S[k].label.toLowerCase()}, not by the ${ref(S[k]).label.toLowerCase()}`, gcol(S[k]));
        item("shared", "run by both");
        item("off", "not run");
      }
    }

    buttons.forEach((b, k) => {
      b.addEventListener("mouseenter", () => { preview = k; show(); });
      b.addEventListener("mouseleave", () => { preview = null; show(); });
      b.addEventListener("focus", () => { preview = k; show(); });
      b.addEventListener("blur", () => { preview = null; show(); });
      b.addEventListener("click", () => { selected = selected === k ? null : k; preview = null; show(); });
    });

    // the tooltip of a function, and the settings that run it
    right.addEventListener("mouseover", ev => {
      const dot = ev.target.closest(".cm-dot");
      if (!dot) return;
      const f = F[+dot.dataset.i];
      tip.textContent = "";
      tip.append(el("code", "", f.name), el("span", "cm-tip-loc", ` src/${f.file}:${f.line}`));
      if (f.doc) tip.append(el("div", "cm-tip-doc", f.doc));
      tip.append(el("div", "cm-tip-used", f.used.length === NS ? "run by every setting" : f.used.length === 0 ?
        "run by none of these settings" : "run by: " + f.used.map(id => S[sidx[id]].label).join(", ")));
      tip.style.display = "block";
      const r = dot.getBoundingClientRect(), w = tip.offsetWidth;
      tip.style.left = Math.max(8, Math.min(r.left + r.width / 2 - w / 2, window.innerWidth - w - 8)) + "px";
      tip.style.top = (r.bottom + 6) + "px";
      f.used.forEach(id => buttons[sidx[id]].classList.add("uses"));
    });
    right.addEventListener("mouseout", ev => {
      if (!ev.target.closest(".cm-dot")) return;
      tip.style.display = "none";
      buttons.forEach(b => b.classList.remove("uses"));
    });
    window.addEventListener("scroll", () => { tip.style.display = "none"; }, { passive: true });

    buildEdges();
    show();
  }

  // the data file (data-src) sits beside this script, whatever the depth of the page
  const base = document.currentScript.src.replace(/[^/]*$/, "");
  function init() {
    document.querySelectorAll("div.codemap[data-var]").forEach(el => {
      const go = () => { const d = window[el.dataset.var]; if (d) setup(el, d); };
      if (window[el.dataset.var] || !el.dataset.src) return go();
      const s = document.createElement("script"); s.src = base + el.dataset.src; s.onload = go;
      document.head.append(s);
    });
  }
  document.readyState === "loading" ? document.addEventListener("DOMContentLoaded", init) : init();
})();
