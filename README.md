# Defect Descriptor

<p align="center">
  <strong>Configurational-force and mixed-mode fracture descriptors directly from full-field experimental maps.</strong>
</p>

<p align="center">
  <a href="https://github.com/Shi2oon/Defect_Descriptor/blob/main/LICENSE"><img alt="License: MIT" src="https://img.shields.io/badge/License-MIT-green.svg"></a>
  <img alt="MATLAB" src="https://img.shields.io/badge/MATLAB-100%25-orange.svg">
  <a href="https://doi.org/10.1007/s00366-025-02262-5"><img alt="Paper DOI" src="https://img.shields.io/badge/DOI-10.1007%2Fs00366--025--02262--5-blue.svg"></a>
</p>

**Defect Descriptor** is a MATLAB toolbox for extracting fracture and defect-field descriptors from measured or synthetic full-field maps. It calculates configurational-force quantities and mixed-mode stress intensity factors from local strain, displacement-gradient, deformation-gradient, or DIC displacement data.

The code is designed for cases where the defect geometry, boundary conditions, or applied load are not easily reduced to a standard fracture-mechanics specimen. Instead of starting from a predefined finite-element model, the workflow starts from measured fields around the defect.

## What the toolbox calculates

| Quantity | Meaning | Typical use |
|---|---|---|
| `$J$` | Energy-release/configurational-force measure | Crack-driving force and contour convergence checks |
| `$M$` | Configurational-force descriptor | Directional defect severity and horizontal/vertical force components |
| `$K_I$` | Mode-I opening stress intensity factor | Tensile crack opening |
| `$K_{II}$` | Mode-II sliding stress intensity factor | In-plane shear |
| `$K_{III}$` | Mode-III tearing stress intensity factor | Anti-plane shear, when the input field supports it |
| Direction sweep | `$J$`, `$M$`, `$K_I$`, `$K_{II}$`, `$K_{III}$` over trial crack directions | Finding the energetically preferred crack or defect direction |

The main solver is:

```matlab
[K, KI, KII, KIII, J, M, Maps] = M_J_KIII_2D(Data, Prop);
```

## Why use it?

Classical fracture analysis normally needs a clean geometry, known loading, and carefully defined contours. Experimental maps from DIC, stereo-DIC, HR-EBSD, or non-standard micro-mechanical tests are rarely that tidy.

This toolbox is useful when you want to:

- extract `$J$`, `$M$`, `$K_I$`, `$K_{II}$`, and `$K_{III}$` directly from local field data;
- analyse non-standard cracks, tortuous cracks, and microstructural defects;
- work with DIC, stereo-DIC, deformation-gradient, displacement-gradient, or HR-EBSD-derived fields;
- perform a direction sweep to identify the maximum-energy defect direction;
- validate the formulation using synthetic calibration fields with known input `$K_I$`, `$K_{II}$`, and `$K_{III}$`.

## Repository layout

```text
Defect_Descriptor/
├── Data/                         # Example datasets
├── functions/                    # Core MATLAB routines
├── input_desk_Validation.m       # Synthetic-field validation examples
├── input_desk_DIC.m              # 2D DIC, stereo-DIC, and elastoplastic examples
├── input_desk_xEBSD.m            # HR-EBSD/xEBSD example workflow
├── input_Direction_Sweep.m       # Direction sweep from -90° to 90°
├── LICENSE                       # MIT licence
└── README.md
```

## Quick start

Clone the repository and run one of the input desks from the repository root:

```matlab
clc; clear; close all
addpath(genpath([pwd '\functions']));
```

Then choose the workflow that matches your input data.

## Example 1: validate using synthetic fields

Use this first. It is the safest way to confirm that the code path is working before using experimental data.

```matlab
clc; clear; close all
addpath([pwd '\functions'])

% Generate a synthetic field with known mixed-mode input
[~, ~, alldata, Prop] = Calibration_2DKIII(3, 1, 2);

% U means displacement-gradient input
Prop.Operation = 'U';

[K, KI, KII, KIII, J, M, Maps] = M_J_KIII_2D(alldata, Prop);
```

You can also validate the deformation-gradient and DIC-displacement input formats through `input_desk_Validation.m`.

## Example 2: DIC or stereo-DIC data

Use `input_desk_DIC.m` for measured displacement fields.

```matlab
clc; clear; close all
addpath(genpath([pwd '\functions']));

Prop.E = 210e9;                  % Young's modulus [Pa]
Prop.nu = 0.30;                  % Poisson's ratio [-]
Prop.units.St = 'Pa';            % Stress unit
Prop.units.xy = 'mm';            % Coordinate unit: 'm', 'mm', 'um', or 'nm'
Prop.stressstat = 'plane_stress';% 'plane_stress' or 'plane_strain'
Prop.Operation = 'DIC';          % Raw DIC displacement data

DataDirect = fullfile(pwd, 'Data', '1KI-2KII-3KII_Data.dat');
Data = importdata(DataDirect);

[K, KI, KII, KIII, J, M, Maps] = M_J_KIII_2D(Data.data, Prop);
```

### Elastoplastic example

For Ramberg-Osgood type behaviour, define the additional material parameters before calling the solver:

```matlab
Prop.E = 210e9;
Prop.nu = 0.30;
Prop.yield = 4e9;          % Yield stress [Pa]
Prop.Yield_offset = 1.24;  % Hardening coefficient
Prop.Exponent = 26.67;     % Hardening exponent
Prop.units.St = 'Pa';
Prop.units.xy = 'm';
Prop.stressstat = 'plane_stress';
Prop.Operation = 'DIC';
```

Check the units carefully. Most bad results from this workflow come from mixing metres, millimetres, micrometres, nanometres, Pa, and MPa.

## Example 3: HR-EBSD / xEBSD data

Use `input_desk_xEBSD.m` when the strain or displacement-gradient field has been prepared through an xEBSD-style workflow.

```matlab
clc; clear; close all
addpath(genpath([pwd '\functions']))

filename = [pwd '\Data\Crack_in_Si_XEBSD'];

[Maps, alldata] = GetGrainData(filename);
M_J_KIII_2D(alldata, Maps);
```

For CrossCourt or other HR-EBSD pipelines, prepare the data into the same structure expected by the other input desks. The example is useful, but it should not be treated as a universal HR-EBSD importer.

## Example 4: direction sweep

`input_Direction_Sweep.m` rotates the local field from `-90°` to `90°` and recalculates the descriptors at each angle.

```matlab
crack_angles = -90:1:90;

for i = 1:length(crack_angles)
    theta = crack_angles(i);

    R = [cosd(theta) -sind(theta) 0;
         sind(theta)  cosd(theta) 0;
         0            0           1];

    % Rotate the tensor field, then solve
    [K, KI, KII, KIII, J, M] = M_J_KIII_2D(RotatedAlldata, Prop);
end
```

The output table can include:

```text
Theta (deg)
J (J/m2)
KI (MPa m^-0.5)
KII (MPa m^-0.5)
KIII (MPa m^-0.5)
J1 (J/m2)
J2 (J/m2)
M1 (J/m)
M2 (J/m)
```

Use this workflow when the crack direction is uncertain, when the defect is not a straight crack, or when the most energetically favourable direction is part of the analysis.

## Input requirements

This part is not optional. The method is sensitive to field formatting.

### Grid requirements

Your data should be:

- square: the number of data points in `$x$` should equal the number of data points in `$y$`;
- uniformly spaced: the spacing in `$x$` and `$y$` should be equal;
- centred: the crack or defect should be placed at the centre of the map;
- consistently scaled: coordinates, stresses, strains, and material constants must use compatible units.

If the map is not square, not uniform, or badly centred, do not expect physically meaningful `$J$`, `$M$`, or `$K$` values.

### Accepted operation modes

| `Prop.Operation` | Input type | Description |
|---|---|---|
| `'DIC'` | Displacement map | Raw 2D or stereo-DIC displacement field |
| `'U'` | Displacement-gradient map | Gradient-based input, useful for validation or HR-EBSD-type data |
| `'F'` | Deformation-gradient map | Deformation-gradient input |

### Strain-map vector format

One accepted format is a vectorised map with coordinates followed by tensor components:

```matlab
Maps = [X(:)  Y(:)  Z(:)  ...
        E11(:) E12(:) E13(:) ...
        E21(:) E22(:) E23(:) ...
        E31(:) E32(:) E33(:)];
```

For a strictly 2D strain map, set the unavailable out-of-plane components to zero. Do not leave missing components as undefined values.

## Material properties

For isotropic elasticity:

```matlab
Prop.E  = 210e9;          % Young's modulus [Pa]
Prop.nu = 0.30;           % Poisson's ratio [-]
Prop.stressstat = 'plane_stress'; % or 'plane_strain'
```

For anisotropic elasticity, use a stiffness matrix instead of `E` and `nu` when the relevant code path requires it:

```matlab
Prop.Stiffness = C;       % stiffness matrix, in Pa
```

For elastoplastic examples:

```matlab
Prop.yield = 4e9;
Prop.Yield_offset = 1.24;
Prop.Exponent = 26.67;
```

Be explicit about the stress state. A plane-stress result and a plane-strain result are not interchangeable.

## Outputs

The solver returns descriptor structures rather than a single scalar:

```matlab
[K, KI, KII, KIII, J, M, Maps] = M_J_KIII_2D(Data, Prop);
```

Typical fields include:

```matlab
J.true        % representative J value
J.div         % scatter/convergence measure
J.Raw         % raw contour values
J.direction_true
J.maxJ_true

KI.true       % representative K_I
KII.true      % representative K_II
KIII.true     % representative K_III

M.true        % M-integral components
M.div         % M-integral scatter/convergence measure
```

Exact structure fields can vary between workflows. Inspect `J`, `M`, `KI`, `KII`, and `KIII` after running the example closest to your data.

## Recommended workflow

```mermaid
flowchart LR
    A[Full-field data] --> B{Input type}
    B -->|DIC or stereo-DIC| C[input_desk_DIC.m]
    B -->|Synthetic validation| D[input_desk_Validation.m]
    B -->|HR-EBSD or xEBSD| E[input_desk_xEBSD.m]
    B -->|Unknown defect direction| F[input_Direction_Sweep.m]

    C --> G[Check units, grid spacing, centring]
    D --> G
    E --> G
    F --> G

    G --> H[M_J_KIII_2D]
    H --> I[$J$, $M$, $K_I$, $K_{II}$, $K_{III}$]
    I --> J[Contour convergence and direction analysis]
```

## Common problems

| Problem | Likely cause | Fix |
|---|---|---|
| Results are unstable across contours | Crack not centred, noisy data, poor contour region | Re-centre the map, crop the field, smooth only where justified, check contour convergence |
| `$K_I$`, `$K_{II}$`, or `$K_{III}$` has the wrong order of magnitude | Unit mismatch | Check `Prop.units.xy`, `Prop.units.St`, `Prop.E`, and data units |
| Solver fails or gives meaningless values | Non-square or non-uniform grid | Resample to a square, uniformly spaced map |
| HR-EBSD example does not run | Missing xEBSD helper function | Add the required `GetGrainData` function or prepare the data manually |
| Direction sweep is slow | One-degree sweep over a large field | Start with `-90:5:90`, then refine around the peak direction |

## Good practice before publication

Before reporting values from experimental data:

1. Run `input_desk_Validation.m` and confirm that the synthetic recovery behaves as expected.
2. Plot the input field and verify the crack or defect centre.
3. Check unit consistency line by line.
4. Report the stress state: plane stress, plane strain, or anisotropic stiffness treatment.
5. Report contour convergence, not only the final scalar value.
6. For direction sweeps, report the angular step and the criterion used to select the preferred direction.

## Citation

If you use this code, cite the associated paper:

```bibtex
@article{Koko2026DefectDescriptor,
  title   = {Bridging experiments and defects' mechanics: a data-driven toolbox for configurational force analysis},
  author  = {Koko, A. and co-authors},
  journal = {Engineering with Computers},
  year    = {2026},
  doi     = {10.1007/s00366-025-02262-5}
}
```

Please check the final author list and bibliographic details from the publisher page before submitting a paper.

## Related project

For DIC-to-Abaqus workflows, see:

- [DIC2ABAQUS](https://github.com/Shi2oon/DIC2ABAQUS)

`DIC2ABAQUS` reconstructs DIC displacement fields inside Abaqus for fracture-mechanics analysis. `Defect Descriptor` instead focuses on extracting defect descriptors directly from full-field maps.

## Licence

This repository is released under the MIT Licence. See [`LICENSE`](LICENSE).
