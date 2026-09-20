# Turning `CrystalVibrationAppDesigner.m` into a `.mlapp`

`CrystalVibrationAppDesigner.m` is `CrystalVibrationApp.m` rewritten the way App Designer
writes its own code: `matlab.apps.AppBase` subclass, public component properties, callbacks
with `(app, event)` signatures wired by `createCallbackFcn`, a generated `createComponents`,
and the standard constructor/`delete` pair.

A `.mlapp` is a binary (zipped) file, so it cannot be written out by hand — it has to be
created by App Designer. This file is the source you copy into it.

## Option A — you may not need a `.mlapp` at all

```matlab
CrystalVibrationAppDesigner
```

The class runs exactly like `CrystalVibrationApp` does. A `.mlapp` only buys you the
Design View editor; it is not required to run, share, or package an app.

## Option B — build the `.mlapp`

### 1. Create the app

App Designer -> **New -> Blank App**, then **Save As** `CrystalVibrationAppDesigner.mlapp`
in this folder (`Version-8`). The class name follows the file name, so if you save it under
a different name, rename the constructor in Code View to match.

### 2. Lay out the components (Design View)

Drag the tree below in, top to bottom, and set each **name** in the Component Browser
exactly as written. The names matter: App Designer derives the callback names from them,
so with these names its generated stubs line up with the code you are pasting.

```
UIFigure                       Position [60 60 1500 880], Name 'Crystal Vibration Explorer', Color white
└─ MainGrid            GridLayout   ColumnWidth {340,'1x'}  RowHeight {'1x',24}  RowSpacing 6  ColumnSpacing 8
   ├─ ControlPanel     Panel        row 1, col 1   Title 'Model & Animation'
   │  └─ ControlGrid   GridLayout   ColumnWidth {'1x'}  Scrollable on  Padding [10 10 10 10]  RowSpacing 6
   │                                RowHeight {24,32,22,24,22,24,22,24,38,26,22,24,22,50,22,50,34,'1x'}
   │     ├─ CrystalHeaderLabel     Label          row 1    'Crystal & bond network' (bold)
   │     ├─ LoadFileButton         Button         row 2    'Load Crystal File...'
   │     ├─ SupercellLabel         Label          row 3
   │     ├─ NEditField             NumericEditField row 4  Limits [1 1], Editable off, Value 1
   │     ├─ kLabel                 Label          row 5
   │     ├─ kEditField             NumericEditField row 6  Limits [0.001 Inf], Value 1
   │     ├─ BondNetworkLabel       Label          row 7
   │     ├─ cutoffDropDown         DropDown       row 8    3 items, default 'Auto (grows shells until rigid)'
   │     ├─ BuildButton            Button         row 9    'Build & Solve'
   │     ├─ AnimationHeaderLabel   Label          row 10   'Mode animation' (bold)
   │     ├─ ModeLabel              Label          row 11
   │     ├─ ModeDropDown           DropDown       row 12   Enable off
   │     ├─ AmplitudeLabel         Label          row 13
   │     ├─ AmplitudeSlider        Slider         row 14   Limits [0 1], Value 0.25
   │     ├─ SpeedLabel             Label          row 15
   │     ├─ SpeedSlider            Slider         row 16   Limits [0.1 10], Value 2
   │     └─ TransportGrid          GridLayout     row 17   ColumnWidth {'1x','1x'}, Padding 0
   │        ├─ PlayPauseButton     Button         col 1    'Play', Enable off
   │        └─ StopButton          Button         col 2    'Stop', Enable off
   ├─ TabGroup         TabGroup     row 1, col 2
   │  ├─ LatticeTab              Tab   'Lattice & Modes'
   │  │  └─ LatticeGrid          GridLayout   ColumnWidth {'5x','4x'}  RowHeight {'1x'}
   │  │     ├─ PlotAx            UIAxes  row 1, col 1
   │  │     └─ SpecAx            UIAxes  row 1, col 2
   │  ├─ ShearTab                Tab   'Shear Analysis'
   │  │  └─ ShearGrid            GridLayout   ColumnWidth {340,'1x','1x'}  RowHeight {'1x','1x'}
   │  │     ├─ ShearControlPanel Panel   rows 1-2, col 1   Title 'Shear Control Panel'
   │  │     │  └─ ShearControlGrid  GridLayout  ColumnWidth {'1x'}  Scrollable on
   │  │     │                       RowHeight {22,50,34,44,26,28,34,34,34,34,96,'1x'}
   │  │     │     ├─ ShearStrainLabel    Label   row 1
   │  │     │     ├─ ShearSlider         Slider  row 2   Limits [0 0.5], Value 0
   │  │     │     ├─ RelaxButton         Button  row 3
   │  │     │     ├─ RelaxStatusLabel    Label   row 4   WordWrap on, VerticalAlignment top
   │  │     │     ├─ V7HeaderLabel       Label   row 5   (bold)
   │  │     │     ├─ GammaRowGrid        GridLayout row 6  ColumnWidth {'1x',110}, Padding 0
   │  │     │     │  ├─ GammaMaxLabel    Label             col 1
   │  │     │     │  └─ GammaMaxEdit     NumericEditField  col 2  Limits [0.05 1.5], Value 0.5
   │  │     │     ├─ PhononPlotsButton   Button  row 7
   │  │     │     ├─ BondPlotsButton     Button  row 8
   │  │     │     ├─ DispersionButton    Button  row 9
   │  │     │     ├─ AlphaRowGrid        GridLayout row 10  ColumnWidth {'1x',80}, Padding 0
   │  │     │     │  ├─ AlphaLabel       Label             col 1  WordWrap on
   │  │     │     │  └─ AlphaEdit        NumericEditField  col 2  Limits [1 20], Value 4
   │  │     │     └─ AnalysisStatusLabel Label   row 11  WordWrap on, VerticalAlignment top
   │  │     ├─ UIAxes3D_NoForce      UIAxes  row 1, col 2
   │  │     ├─ UIAxes3D_WithForce    UIAxes  row 1, col 3
   │  │     ├─ UIAxesMode_NoForce    UIAxes  row 2, col 2
   │  │     └─ UIAxesMode_WithForce  UIAxes  row 2, col 3
   │  ├─ PhononTab               Tab   'Phonon Softening'     (leave empty)
   │  ├─ BondTab                 Tab   'Bond Breaking'        (leave empty)
   │  └─ DispersionTab           Tab   'Dispersion Surface'   (leave empty)
   └─ StatusLabel      Label        row 2, cols 1-2
```

The last three tabs stay empty on purpose: `plotShearSoftening`, `plotShearBonds` and
`brillouinShearSurface` build their contents into them at run time.

For the exact colours, tooltips and default text of every component, read `createComponents`
in `CrystalVibrationAppDesigner.m` — it is written in App Designer's own generated style, so
it reads as a line-by-line checklist. You cannot paste it into App Designer: App Designer
generates that section itself and keeps it read-only in Code View.

### 3. Paste the private properties

Code View -> **Property -> Private Property**. App Designer makes an empty
`properties (Access = private)` block; replace its contents with the block of the same name
from `CrystalVibrationAppDesigner.m` (`latticeVectors` through `AnimFrameRateHz`).

### 4. Paste the helper functions

Code View -> **Function -> Private Function**, once per helper. Rename the stub, then paste
the body. The helpers are, in file order:

`resetAnalysisTabs`, `clearAnalysisTab`, `setShearControlsEnabled`, `refreshShearViewIfNeeded`,
`runAnalysisPlots`, `updateShearAnalysis`, `populateModeDropdown`, `getSelectedModeIndex`,
`plotFrequencySpectrum`, `plotStaticLattice`, `highlightModeOnSpectrum`, `speedToPeriod`,
`animationTick`, `stopTimerIfRunning`, `setStatus`, `showAlert`, `makeHoverLabels`.

Two of them have non-default signatures — keep them as written:
`clearAnalysisTab(~, tab, kind)` and `labels = makeHoverLabels(~, elementNames, masses)`.

### 5. Add the callbacks

Right-click the component -> **Callbacks -> Add ... callback**. With the names from step 2,
App Designer generates precisely these stubs; paste the matching body into each.

| Component | Callback to add | Generated name |
| --- | --- | --- |
| UIFigure | CloseRequestFcn | `UIFigureCloseRequest` |
| UIFigure | StartupFcn | `startupFcn` |
| LoadFileButton | ButtonPushedFcn | `LoadFileButtonPushed` |
| BuildButton | ButtonPushedFcn | `BuildButtonPushed` |
| ModeDropDown | ValueChangedFcn | `ModeDropDownValueChanged` |
| PlayPauseButton | ButtonPushedFcn | `PlayPauseButtonPushed` |
| StopButton | ButtonPushedFcn | `StopButtonPushed` |
| TabGroup | SelectionChangedFcn | `TabGroupSelectionChanged` |
| ShearSlider | ValueChangedFcn | `ShearSliderValueChanged` |
| RelaxButton | ButtonPushedFcn | `RelaxButtonPushed` |
| PhononPlotsButton | ButtonPushedFcn | `PhononPlotsButtonPushed` |
| BondPlotsButton | ButtonPushedFcn | `BondPlotsButtonPushed` |
| DispersionButton | ButtonPushedFcn | `DispersionButtonPushed` |
| AlphaEdit | ValueChangedFcn | `AlphaEditValueChanged` |

`startupFcn` is added from the **app** (right-click `app.UIFigure` in the Component Browser),
not from a control.

### 6. Run it

`.mlapp` must sit in `Version-8` (or have it on the path): the app calls `autoRigidShells`,
`buildDynamicalMatrix`, `jacobiEigenSolver`, `applyShearForce`, `relaxShearLattice`,
`plotShearSoftening`, `plotShearBonds`, `brillouinShearSurface`, `shearPlotCanvas`,
`plotLattice`, `plotLatticeWithGradient`, `atomHover`, `getElementDetails`, `parseCIFFile`
and `parseCrystalFile`.

## Differences from `CrystalVibrationApp.m`

The GUI and the numerics are identical. Only the plumbing changed, to match what App
Designer expects:

- `createComponents` is flat rather than split into `buildControlColumn` / `buildLatticeTab` /
  `buildShearTab` / `buildShearControls`, because App Designer generates one flat function.
  Every container therefore needs a name, so the layout grids and the static labels are now
  named properties (`MainGrid`, `ControlGrid`, `kLabel`, ...) instead of local variables.
- Callbacks take `(app, event)` and are wired with `createCallbackFcn` instead of anonymous
  functions. `updateShearAnalysis` now reads `app.ShearSlider.Value` through
  `ShearSliderValueChanged` rather than taking the value from an event struct, and the
  inline `AlphaEdit` handler became the named `AlphaEditValueChanged`.
- `TabSelectionChanged` is renamed `TabGroupSelectionChanged`, and `openDispersionSurface`
  is renamed `DispersionButtonPushed`, so both match App Designer's naming.
- `makeHoverLabels` moved from a local function after the `classdef` to a private method
  (App Designer has no place for local functions), so calls are now `app.makeHoverLabels(...)`.
- `onClose` became the `UIFigureCloseRequest` callback, and `AnimFrameRateHz` moved from a
  `Constant` block into the ordinary private properties, since App Designer's property
  editor only creates plain public/private blocks.

One consequence worth knowing: App Designer generates `delete(app)` and keeps it read-only,
so the animation timer is stopped in `UIFigureCloseRequest` (closing the window) rather than
in `delete`. If you tear the app down with `delete(app)` from the command line instead of
closing the window, stop the timer yourself first.
