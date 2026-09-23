# The App Designer app: `app1.mlapp`

`app1.mlapp` is the Crystal Vibration Explorer built in App Designer (R2024b). It uses
the same helper functions as `CrystalVibrationApp.m`, the version of the app written
entirely in code. It also has the same numerics and the same five tabs: Lattice and
Modes, Shear Analysis, Phonon Softening, Bond Breaking, Dispersion Surface.

## Running it

From the `Version-8` folder:

```matlab
launch_GUI          % opens app1
app1                % same thing
CrystalVibrationApp % the pure-code version, still works
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

## Component names

Most components use the names from `CrystalVibrationAppDesigner.m`. The labels, panels and
top-level grid use the names App Designer gave them automatically. The code only uses one of
these renamed components, `MorseAlphaLabel`.

| `CrystalVibrationAppDesigner.m` | `app1.mlapp` |
| --- | --- |
| `MainGrid` | `GridLayout` |
| `ControlPanel` | `ModelandAnimationPanel` |
| `SupercellLabel` | `SuperCellEditFieldLabel` |
| `kLabel` | `KEditFieldLabel` |
| `BondNetworkLabel` | `BondNetworkDropDownLabel` |
| `AnimationHeaderLabel` | `ModeAnimationLabel` |
| `ModeLabel` | `VibrationalModeDropDownLabel` |
| `AmplitudeLabel` | `AmplitudeSliderLabel` |
| `SpeedLabel` | `PlaybackSpeedSliderLabel` |
| `ShearStrainLabel` | `EngineeringShearStrainGammaLabel` |
| `GammaMaxLabel` | `MaxmimumSheerStrainGammaLabel` |
| `AlphaLabel` | `MorseAlphaLabel` |

`app1.mlapp` also has two empty labels, `Label` and `Label_2`, that the code never uses.

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
