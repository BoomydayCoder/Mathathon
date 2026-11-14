# Local verification of the pinned Romik Lean project

**Status: local compilation and exact theorem/axiom checks passed on
2026-11-14 UTC.** The full pinned project built successfully (8995 jobs, exit
status 0). Both the repository's statement check and the independent fully
qualified statement check compiled with exit status 0. The separate project
kernel replay is still running; no completed replay result is claimed yet.

Source: `https://github.com/danromik/zeta5-irrationality`, commit
`fc2fefa77076e22e7e80b4c101ded5a98d53e14d`, Git tree
`738083ba7638325055bb2a37b795c1a57d463327`.

Toolchain: Lean 4.34.0, compiler commit
`293d5d0c0c3f3dded4688b3ccd6a33939ac5102b`; Mathlib commit
`5ed2965256430c3649e86755f9576b54eca72435`, as pinned by the repository manifest.
The complete manifest and binary hashes are recorded in
`../verification/romik-lean/metadata.json`.

All downloaded source, toolchains and build artifacts are under the ignored
`codex/.cache/lean-romik/`. Only `.cache/` was added to `codex/.gitignore`.
No source file in the pinned project has been edited.

The full source has 70 files under `Zeta5/`, totaling 28,976 lines. The four
files previously archived in `codex/sources/` match byte-for-byte. Source scans
find the explicit PNT axiom and no use of `sorry`, `admit`, `native_decide`,
`implemented_by`, `extern` or unsafe declarations in the project proof code;
mentions of such terms in comments do not count as declarations. These scans supplement the completed compiled axiom report; they are not its
replacement.

The declared external axiom is the usual prime number theorem:

```lean
Tendsto (fun x : ℝ => Chebyshev.theta x / x) atTop (nhds 1)
```

The requested final statement is irrationality of the ordinary real series
\(\sum_{v\ge0}(v+1)^{-5}\). An additional check using fully qualified names
was compiled successfully from `../verification/romik-lean/Check.lean`.

Reproduction commands, from the pinned checkout, use the following environment:

```sh
export ELAN_HOME=/workspace/Mathathon/codex/.cache/lean-romik/elan
export MATHLIB_CACHE_DIR=/workspace/Mathathon/codex/.cache/lean-romik/mathlib-cache
export PATH="$ELAN_HOME/bin:$PATH"
export LEAN_NUM_THREADS=2
lake exe cache get
lake --no-cache build
lake env lean cert/Check.lean
lake env lean /workspace/Mathathon/codex/verification/romik-lean/Check.lean
```

All displayed commands completed successfully. `--no-cache` on the build
prevented downloading project build artifacts; the already fetched Mathlib
objects remained available as dependencies.

The installed `LeanChecker.lean` source was inspected. Default
`lake env leanchecker --verbose Zeta5` replays each project module's declarations
against its imported environment. It does not independently recheck every
Mathlib declaration. The separate `--fresh Zeta5` mode replays imported constants
as well into an empty environment. Both use Lean's own kernel, not a second
independent kernel implementation. The first replay, started at 09:54:25 UTC, lost its wrapper/result capture.
Its old process was later a zombie with Linux exit-code field zero and no
recorded OOM, but the log contained only the initial banner and no result JSON.
The cause of the missing wrapper is not established; this incomplete record is
not used as a completed replay certificate. The evidence is preserved in
`kernel-replay-attempt1.log` and `kernel-replay-attempt1-status.json`.

A second replay started at 10:17:08 UTC with one Lean worker, PTY-backed logging,
a periodically updated result JSON, and a 45-minute bound. It is running while
the reviewing agent remains active. No completed replay result is claimed yet.


## Actual compiled statement and axiom output

The independent check compiled this declaration:

```lean
example : _root_.Irrational
    (_root_.tsum (fun v : _root_.Nat =>
      (1 : _root_.Real) / ((v : _root_.Real) + 1) ^ 5)) :=
  Zeta5.zeta5_irrational
```

The actual output was:

```text
def Zeta5.zeta5 : ℝ :=
∑' (v : ℕ), 1 / (↑v + 1) ^ 5
axiom Zeta5.Axioms.chebyshev_theta_asymptotic :
  Filter.Tendsto (fun x => Chebyshev.theta x / x) Filter.atTop (nhds 1)
'Zeta5.zeta5_irrational' depends on axioms: [propext,
 Classical.choice,
 Quot.sound,
 Zeta5.Axioms.chebyshev_theta_asymptotic]
```

No `sorryAx` or additional project axiom occurs. The 17 build warnings were
unused-variable/unused-section-variable linters, matching the categories stated
by the author; there was no `sorry` warning. The checkout remains clean and
unchanged at the pinned commit. The pinned commit differs from the author's
previously reported tested commit `a87234f` only in `docs/CERTIFICATION.md`;
the current successful local build does not rely on that historical report.

Saved evidence:

- [Full local build log](../verification/romik-lean/build.log)
- [Repository statement check](../verification/romik-lean/repository-check.log)
- [Independent statement/axiom check](../verification/romik-lean/independent-check.log)
- [Exact source, compiler and log hashes](../verification/romik-lean/metadata.json)
- [Kernel replay log, currently incomplete](../verification/romik-lean/kernel-replay.log)

The mathematical external input is PNT, an established theorem, but its proof is
not formalized within this pinned development. The locally checked formal claim
is therefore the implication from that PNT statement to irrationality of the
ordinary real ζ(5). Lean's standard logical axioms and kernel remain trusted;
Mathlib dependencies were obtained from the pinned compiled cache.

This check concerns Romik's unmodified formalization of Fauzan's argument. It
does not make the new equilibrium construction or Python certificates in this
research branch Lean-verified. Those have their own mathematical and numerical
reviews. Nor is local compilation evidence of a historical novelty claim or
external scholarly acceptance.
