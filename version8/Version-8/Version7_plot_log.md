# Version 7 – how every plot is built (what comes from your code, what is assumed)

Run from the app (Shear Analysis window → two green/orange buttons) or from `runPhononShearAnalysis.m`
(set `cifFile`, saves PNGs to `shear_figures/`). Both call the same functions, so results are identical.

## 0. Pipeline

| Step | File | Origin |
|---|---|---|
| Read crystal | `parseCIFFile.m` / `getBuiltinCrystalDef.m` | **your code** |
| Atoms, bonds, masses | `autoRigidShells.m` -> `generateLatticeGeneral.m` | **your code** (bonds = `[i j d0 ux uy uz]`, one spring per neighbour image, same k for every bond) |
| gamma = 0 spectrum | `buildDynamicalMatrix.m` (checked, see §4) | **your code** |
| Shear geometry | `phononShearSweep.m` | **new**, but uses your convention u_x = gamma*z (same as `applyShearForce.m`) |
| 3 model levels | `phononShearSuite.m` | **new** |
| Plots | `plotShearSoftening.m` (P1-P6), `plotShearBonds.m` (Q1-Q2) | **new** |
| Transport proxies | `phononTransportMetrics.m` | **new** (assumptions in §3) |
| Helpers | `assembleShearHessian.m`, `shearBondModel.m` | **new** |

Existing functions reused inside the new code: `extractModeShape` (P4), `centralDiffDerivative` (P5, P6),
`jacobiEigenSolver` (optional flag `useJacobi`, used only for the cross-check). `applyShearForce`, `computeShearForces`,
`buildDynamicalMatrix` are used in `runPhononShearAnalysis.m` as validation references.

## 1. The three model levels (what each curve means)

1. **harmonic, affine** – identical to your sliders: every bond keeps stiffness k, only its *direction* rotates with the shear
   (u_x = gamma*z applied to every atom). No relaxation, no breaking. Verified equal to `applyShearForce.m` to 1e-15.
2. **harmonic + pre-stress + relaxation** – adds (a) the geometric stiffness of stretched/compressed bonds, the exact second
   derivative of the bond energy: block = V'' u u' + (V'/d)(I - u u'), (b) internal (non-affine) atom relaxation.
3. **Morse + pre-stress + relaxation + rupture** (main model) – bonds soften as they stretch and are removed when they
   pass their force maximum.

## 2. Data → figure

Common data per strain step (all computed by `phononShearSweep.m` from **your** pos/bonds/masses/k):
bond vectors r = r0 + gamma*r0_z*x-hat (+ relaxation w_j - w_i), bond stretch (d-d0)/d0, Hessian K -> D = M^-1/2 K M^-1/2,
translations projected out (3 exact zero modes), eigenvalues lambda = omega^2, lowest 6 eigenvectors, virial shear stress tau_xz.

| Fig | What is drawn | Data source | Assumed / derived |
|---|---|---|---|
| **P1** softening curve | Left: soft-band RMS(omega) of the lowest 5 % of modes, and mean omega of all modes, vs gamma (main model). Right: soft band for the 3 models. | eigenvalues of D at each gamma (**your** K construction, rotated bond directions) | Soft band = lowest 5 % of modes (min 3): chosen instead of the single lowest mode because that one is fragile for marginally rigid networks (Al2O3 lowest lambda = 8e-4 at gamma=0). Curves stop at the first strain with omega^2<0 (dotted last segment). |
| **P2** DOS shift | Gaussian-broadened DOS at gamma=0 vs last stable gamma; mean omega and soft-band shifts in the title | same eigenvalues | Broadening sigma = 1.5 % of the top frequency (visual only). |
| **P3** spectral map | DOS(omega, gamma) image, lowest mode in red | same eigenvalues | Same broadening; only stable range shown. |
| **P4** softest mode | All atoms projected on xz, size/colour = \|u_i\|, arrows = in-plane displacement, gamma=0 vs last stable gamma | eigenvector of lowest mode -> **`extractModeShape.m`** (un-mass-weighting + normalisation); atom positions = pos + gamma*z + relaxation | Projection onto xz overlaps atoms in 3D crystals (larger atoms drawn behind smaller). |
| **P5** stress vs phonon | tau_xz (left axis); soft-band ratio and sqrt(G_t/G_0) (right axis) | tau = virial sum f*r_x*r_z/(d*Vol) from bond forces (same virial formula as **`computeAtomicStress.m`**, checked against **`computeShearForces.m`**); G_t = d tau/d gamma via **`centralDiffDerivative.m`** on a pchip interpolant | Vol = Natoms x (cell volume / atoms per cell). Shear-wave speed v_s ∝ sqrt(G) is standard elasticity, not computed from k-space. |
| **P6** transport | (a) mean-square displacement, (b) phonon-limited mobility proxies, (c) v_sat and Debye-temperature proxies | eigenvalues (all modes) | see §3 |
| **Q1** bond map | every bond projected on xz at 3 strains (0.5x first-break, at first break, last stable); colour = stretch; black x/dashed = broken | bond list + stretch + `active` flags from the sweep | Breaking rule (§3). |
| **Q2** statistics | (a) % broken bonds vs gamma per bond length; (b) stretch vs n_x n_z at last stable gamma with the first-order line eps = gamma*n_x*n_z | bond unit vectors `bonds(:,4:6)` (**your** list) and computed stretch | Critical-stretch line = ln2/alpha. |

## 3. Everything that is ASSUMED (not from your code)

* **Bond potential beyond Hooke**: Morse, V = D(1-e^{-a(d-d0)})^2, D = k/(2a^2) so stiffness at rest equals your k.
  a = alpha/d0 per bond with **alpha = 4**, i.e. every bond has the same critical stretch ln2/alpha = **17.3 %**.
  Change `alphaMorse` in `phononShearSweep.m` to move it.
* **Bond breaking**: a bond is removed permanently when d - d0 exceeds the Morse force maximum. Applies to *all* springs in your
  network, including the long-range ones `autoRigidShells` adds for rigidity (Al2O3 needs 4 shells), so "bond" here means "spring".
* **Same k for every spring** (your model) – no chemistry-specific force constants.
* **Pre-stress term and non-affine relaxation**: physics added on top of your affine model. Relaxation solves for periodic
  internal displacements w (bond vector = affine + w_j - w_i, exact for wrapped bonds), warm-started from the previous strain,
  Newton with line search; residual force is stored (`relaxResidual`, ~1e-9). It is a different (periodic) relaxation from your
  `relaxShearLattice.m` button, which fixes top/bottom layers.
* **Instability strain** = first gamma with lambda < -1e-8*max(lambda). Beyond it the relaxed structure has left the original
  configuration, so all curves are truncated there.
* **Gamma-point only**: a supercell of N cells samples only the k-points that fold to Gamma. Larger N can reveal softening at
  finite k that N=1 cannot (Al2O3: instability at gamma≈0.31 for N=1, ≈0.17 for N=2).
* **Reduced units**: omega in sqrt(k/m) with the masses from the CIF and the k in the app.
* **Transport proxies (P6)** – scaling laws, not device predictions:
  - <u^2> ∝ (1/N) sum hbar/(2 m omega) coth(hbar omega/2kT) (Debye–Waller);
    the top unstrained mode is mapped to **15.5 THz** to fix hbar*omega/kT; T = 77, 300, 450 K; omega floored at 5 % of its
    unstrained lowest value so unstable modes do not diverge.
  - mobility mu_DW/mu_0 = <u^2>_0/<u^2> (high-T phonon-limited resistivity ∝ <u^2>); mu_ac/mu_0 = G_t/G_0 (deformation-potential
    acoustic mobility ∝ elastic constant).
  - v_sat/v_sat,0 = sqrt(omega_op/omega_op,0) with omega_op = mean of the top 10 % of modes (optical-emission-limited saturation
    velocity); Theta_D ratio from the log-mean frequency.
  - Band-structure effects (which dominate real strained-silicon mobility gains) are **not** included.

## 4. Validation (printed by `runPhononShearAnalysis.m`, Al2O3)

* affine model vs `applyShearForce.m` frequencies: max |d omega| = 7e-16
* affine tau_xz vs `computeShearForces.m` virial: identical (0.047766)
* gamma=0 spectrum vs `buildDynamicalMatrix.m`: 1e-15
* residual force after relaxation: < 1e-9
* the script also runs all four CIFs in `CIF files/` and prints first-break strain, instability strain and soft-band ratio.

## 5. What changed for CIF crystals (why Al2O3 looked erratic before)

1. The single lowest mode of a marginally rigid network (Al2O3 lowest lambda = 8e-4 vs ~0.6 top) jumps around and is stiffened by
   pre-stress ("tensegrity" effect) before it collapses -> replaced by the soft-band RMS.
2. Relaxation used raw Newton steps -> now line-searched with a convergence check.
3. Morse parameter used one absolute length scale, so long springs had unrealistically different critical strains ->
   now the same relative critical stretch for every bond.
4. Points after the instability (relaxed structure on a different branch) were plotted -> now cut off.
5. The y-slice used for the 2D pictures showed only a few atoms for non-cubic cells -> now all atoms/bonds are projected on xz.
