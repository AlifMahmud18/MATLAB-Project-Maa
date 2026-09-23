# The App Designer app: `app1.mlapp`

`app1.mlapp` is the Crystal Vibration Explorer, built in App Designer (R2024b). It has five
tabs: Lattice and Modes, Shear Analysis, Phonon Softening, Bond Breaking and Dispersion
Surface. All the numerics live in the helper `.m` files in this folder.

## Running it

From the `Version-8` folder:

```matlab
launch_GUI   % opens app1
app1         % same thing
```

The app must stay in `Version-8`, or that folder must be on the path. It calls
`autoRigidShells`, `generateLatticeGeneral`, `buildDynamicalMatrix`, `jacobiEigenSolver`,
`applyShearForce`, `relaxShearLattice`, `phononShearSuite`, `plotShearSoftening`,
`plotShearBonds`, `brillouinShearSurface`, `shearPlotCanvas`, `plotLattice`,
`plotLatticeWithGradient`, `atomHover`, `getElementDetails`, `getBuiltinCrystalDef`,
`parseCIFFile` and `parseCrystalFile`.

The class name is `app1` because the file is called `app1.mlapp`. To rename it, use
**Save As** in App Designer, which renames the class too. Then change the last line of
`launch_GUI.m` to match. Renaming the file in a file browser gives an error, because the
class name inside the file no longer matches the file name.

The Phonon Softening, Bond Breaking and Dispersion Surface tabs are empty in Design View on
purpose. `plotShearSoftening`, `plotShearBonds` and `brillouinShearSurface` build their
contents at run time. The three buttons in the Shear Analysis tab trigger them.

## Component names used by the code

The code refers to components by name, so don't rename these in Design View:

- **Controls:** `NEditField`, `kEditField`, `cutoffDropDown`, `ModeDropDown`,
  `AmplitudeSlider`, `SpeedSlider`, `PlayPauseButton`, `StopButton`.
- **Lattice and Modes tab:** `TabGroup`, `LatticeTab`, `PlotAx`, `SpecAx`.
- **Shear Analysis tab:** `ShearTab`, `ShearSlider`, `RelaxButton`, `RelaxStatusLabel`,
  `GammaMaxEdit`, `AlphaEdit`, `MorseAlphaLabel`, `PhononPlotsButton`, `BondPlotsButton`,
  `DispersionButton`, `AnalysisStatusLabel`, `UIAxes3D_NoForce`, `UIAxes3D_WithForce`,
  `UIAxesMode_NoForce`, `UIAxesMode_WithForce`.
- **Plot tabs:** `PhononTab`, `BondTab`, `DispersionTab`.
- **Window and panels:** `UIFigure`, `StatusLabel`, `ModelandAnimationPanel`,
  `ShearControlPanel` (the last two are used only by `applyDarkTheme`).

Everything else, including the two empty labels `Label` and `Label_2`, can be renamed or
deleted freely.

## Standalone scripts

These run without the app:

- **`runPhononShearAnalysis.m`** runs the same shear sweeps as the Phonon Softening and Bond
  Breaking tabs, and saves the plots as PNGs in `shear_figures/`.
- **`runStressStrainSweep.m`** runs the stress-strain analysis. It is the only user of
  `checkPlasticity`, `compositeTrapz`, `computeAtomicStress` and `shearYieldGap`.

## Fixes applied to the code in `app1.mlapp`

The design was kept exactly as drawn: layout, colours and component properties are
unchanged. Three things in the code were fixed:

1. **`AlphaEditValueChanged` crashed.** It wrote to `app.AlphaLabel`, which does not exist
   in this app, so changing Morse alpha raised an error. It now updates `app.MorseAlphaLabel`,
   for example `Morse Alpha (break 17.3%)` for alpha = 4. The label is short enough to fit
   its column. The percentage is the bond stretch at which a Morse bond breaks, ln 2 / alpha.
2. **Closing the window raised an error.** `UIFigureCloseRequest` called
   `app.applyDarkTheme()` *after* `delete(app)`, when the app no longer exists. That call is
   now commented out. `applyDarkTheme` itself is kept, but nothing calls it, because the call
   in `startupFcn` was already commented out. The theme is set in Design View instead.
3. **Status text was unreadable on the dark background.** The code coloured the bottom
   status bar and the Relax status label dark blue, dark green, dark red and grey, the
   colours of the old white layout. They are now light versions of the same colours:
   working = light blue, done = light green, failed = light red, idle = light grey. The grey
   placeholder text in the three empty plot tabs keeps its grey, because those tabs are light.

## Suggested tidy-ups in Design View

These are component properties, so they have to be changed in Design View. Nothing breaks
without them.

- Typos: `Maxmimum Sheer Strain (Gamma)` should be `Maximum Shear Strain (Gamma)`, and
  `3D Dispersion Surfce` should be `3D Dispersion Surface`.
- The Shear Analysis control panel is titled `Panel`, and the window is titled `MATLAB App`.
- The build button says **Solve and Build**, but the status and alert messages say
  "Build & Solve". Rename the button, or leave it; it only affects wording.
- `RelaxStatusLabel` has `WordWrap` off. The residual-force message it shows after
  **Relax** is two sentences long, so it gets cut off. Turn `WordWrap` on.
- `MorseAlphaLabel` starts as plain `Morse Alpha` and shows the break stretch only after
  you change the value. To show it from the start, set its Text to
  `Morse Alpha (break 17.3%)`, which matches the default alpha of 4.
