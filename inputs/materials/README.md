# Material property data

The simulation loads defaults from `<materials_dir>/<Material>/<property>.txt`.
Folders and filenames are case-sensitive. Edit the data files without rebuilding.
Each run loads the files once when its config is read. 
<!-- Files are not watched
while a simulation runs. -->

| Material folder | Default properties |
| --- | --- |
| Graphite | ocv, chemical_potential, exchange_current_density, mobility, conductivity, site_density |
| LFP | ocv, chemical_potential, exchange_current_density, diffusivity, mobility, conductivity, site_density|
| NMC | ocv, chemical_potential, exchange_current_density, diffusivity, conductivity, site_density |
| Carbon | ocv, chemical_potential, exchange_current_density, diffusivity, conductivity, site_density |
| Silicon | ocv, chemical_potential, exchange_current_density, diffusivity, conductivity, site_density |
| Electrolyte | diffusivity |

Unavailable properties are not synthesized. For example, graphite transport
uses mobility, and no graphite diffusivity default is supplied.

## File format and units

A constant file contains one number, optionally with `#` comments:

```text
# Graphite conductivity [S/cm]
3.3
```

A curve contains at least two concentration/value rows:

```text
# stoichiometry OCV [V] -- illustrative only
0.0 0.9
0.5 0.2
1.0 0.1
```

All numbers must be finite. Concentration must be strictly increasing, in [0,1]
for electrodes and nonnegative for electrolyte. Site density requires a scalar.
Diffusivity, mobility, conductivity, and exchange current density must be
nonnegative; site density must be positive. Other property values can be signed.

| Property | Units/convention |
| --- | --- |
| ocv | V |
| exchange_current_density | A/cm² |
| diffusivity | cm²/s |
| conductivity | S/cm |
| site_density | mol/cm³ |
| chemical_potential | Preserves the original solver conventions: NMC/LFP/Silicon use -Frd × OCV; Carbon uses -OCV; Graphite uses its independent legacy data |
| mobility | Original Cahn–Hilliard solver convention; graphite's legacy ×(100 × 2/3) conversion is already included |

<!-- Electrode table concentration is stoichiometry; electrolyte concentration is
mol/cm³. Linear interpolation is used between rows. Default files clamp to the
nearest endpoint outside their range. This changes analytical extrapolation:
all electrode curves now clamp outside [0,1]. The electrolyte table covers
**0–0.01 mol/cm³** (0–10 mol/L); supply a wider table if the run needs it. -->

## Overrides and locating defaults

The makefile embeds the absolute `inputs/materials` directory in the executable.
To relocate the data or choose a different complete default library, use:

```ini
materials_dir = materials
```

That directory is relative to the run config (absolute paths also work).
Alternatively set `MATERIALS_DIR=/path/to/materials` when building. All default
files listed above must exist unless replaced by a property override. An
explicitly selected directory never silently falls back to another library.

Override individual properties using the existing config interface:

```ini
material.Graphite.ocv = table materials/Graphite/ocv.txt clamp
material.Graphite.site_density = table materials/Graphite/site_density.txt
material.Graphite.conductivity = constant 4.0
material.NMC.ocv = table my_measurements/nmc_ocv.txt error
```

Override filenames are arbitrary and relative to the config file. Quote table
paths containing spaces. `clamp` is optional; `error` rejects evaluation outside
the supplied table range. Scalar files have no concentration range.

Default OCV and chemical-potential files are independent: editing one does not
rewrite the other. Keep physically related tables consistent when editing the
library. For compatibility with the earlier override interface, a **config OCV
override** for NMC/LFP/Silicon/Carbon also supplies a scaled chemical potential
unless the config explicitly overrides chemical_potential. Graphite always
keeps its independent chemical-potential data. Diffusivity overrides never
implicitly change mobility.
<!-- 
## Migration provenance and numerical accuracy

The initial files preserve the pre-migration MaterialProperties.cpp models:

- Graphite and LFP retain their original table knots. Graphite exchange current
  includes ×1e-3; graphite mobility includes ×(100 × 2/3). LFP OCV includes the
  `3.4 - original_value` conversion, and chemical potential includes -Frd.
- Carbon exchange current retains the original graphite-derived data, now in
  Carbon's own file. Constant properties are stored without resampling.
- Previously analytic curves were sampled adaptively, checking linear
  interpolation at quarter, half, and three-quarter points. The sampling target
  was 1e-6 relative plus one tenth of the absolute tolerances below.
- Silicon retains the original clipped endpoint behavior at 1e-8 and 1-1e-8,
  with dense sampling near those endpoints. The large endpoint OCV values are
  inherited from the old fit, not new measurements.
- Electrolyte diffusivity includes the current D0=0.00489. Chemical-potential
  tables use Frd=96485.3365. Changing Constants.cpp does not regenerate tables.
- Old source vectors live in Graphite/legacy and LFP/legacy for plotting and
  provenance. The simulation never loads those vectors.

The migration comparison uses 20,000 concentrations plus 12 endpoint/interior
probes for each of 33 properties. Its tolerance is 2e-5 relative plus: OCV and
Carbon chemical potential 1e-7; other chemical potentials 0.01; exchange current
1e-12; conductivity 1e-10; electrode diffusivity 1e-16; electrolyte diffusivity
1e-12; mobility 1e-20; site density 1e-15; LFP chp_value 1e-12. The 396 fixed
pre-migration values are retained in `tests/material_defaults_reference.txt`
and checked by `make test-materials`. These are interpolation checks, not a
claim of full simulation equivalence. -->
