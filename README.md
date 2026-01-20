# ATR Permittivity Extraction

Extract frequency-dependent complex permittivity ε(ω) from Attenuated Total Reflection (ATR) spectroscopy measurements using Drude-Lorentz model fitting.

## Overview

This tool uses inverse modeling to extract the complex permittivity of doped semiconductor thin films from ATR reflectance measurements. It combines:
- **Drude model** for free carrier (plasmonic) contributions
- **Lorentz model** for phonon (lattice vibration) contributions

The permittivity is described by:

```
ε(ω) = ε_∞ - ω_p²/(ω² + iγ_p·ω) + Σ[(S_j·ω_TO,j²)/(ω_TO,j² - ω² - iγ_j·ω)]
```

Where:
- `ε_∞` = high-frequency dielectric constant
- `ω_p` = plasma frequency (free carriers)
- `γ_p` = plasma damping rate
- `S` = oscillator strength (phonons)
- `ω_TO` = transverse optical phonon frequency
- `γ_ph` = phonon damping rate

## Features

- **Automated fitting** using differential evolution global optimization
- **Configurable via JSON** - no code editing required
- **Multiple sample support** - process entire datasets at once
- **Hall effect initialization** - uses carrier concentration and mobility for initial guesses
- **Comprehensive visualization** - fit quality, real and imaginary permittivity
- **Export formats** - text files compatible with TMM, FDTD, and other optical simulation tools

## Installation

### Requirements

```bash
pip install numpy scipy matplotlib
```

Or using conda:

```bash
conda install numpy scipy matplotlib
```

### Clone Repository

```bash
git clone https://github.com/yourusername/atr-permittivity-extraction.git
cd atr-permittivity-extraction
```

## Quick Start

### 1. First Run - Generate Configuration Template

```bash
python extract_permittivity_from_ATR.py
```

This creates a `config.json` template file.

### 2. Edit Configuration

Open `config.json` and fill in your parameters:

```json
{
    "material_constants": {
        "name": "ScN",
        "epsilon_inf": 12.8,
        "omega_TO": 359.0,
        "omega_LO": 677.0
    },
    "atr_setup": {
        "prism_material": "ZnSe",
        "prism_index": 2.4,
        "incident_angle_deg": 60
    },
    "samples": {
        "sample_A": {
            "label": "High-doped p-type ScN",
            "thickness_nm": 370,
            "carrier_concentration_cm3": 1.8e20,
            "mobility_cm2Vs": 10,
            "carrier_type": "p-type",
            "effective_mass_me": 0.8,
            "data_directory": "data/sample_A",
            "data_pattern": "*.dpt",
            "fit_range_cm1": [300, 1000],
            "plot_range_cm1": [200, 1200]
        }
    }
}
```

### 3. Prepare Data Directory

Organize your ATR data files:

```
project/
├── extract_permittivity_from_ATR.py
├── config.json
├── data/
│   ├── sample_A/
│   │   └── measurement.dpt
│   └── sample_B/
│       └── measurement.dpt
└── output/
```

Data files should be two-column format (wavenumber, reflectance):
```
400.5, 0.234
401.0, 0.236
401.5, 0.238
...
```

### 4. Run Extraction

```bash
python extract_permittivity_from_ATR.py
```

Or with a custom config file:

```bash
python extract_permittivity_from_ATR.py my_config.json
```

## Configuration Guide

### Material Constants

These are typically fixed from literature or separate measurements:

```json
"material_constants": {
    "name": "Material name (e.g., ScN, InN, GaN)",
    "epsilon_inf": 12.8,           // High-frequency dielectric constant
    "omega_TO": 359.0,             // TO phonon frequency (cm⁻¹)
    "omega_LO": 677.0              // LO phonon frequency (cm⁻¹) - for plotting only
}
```

### ATR Setup

Your experimental configuration:

```json
"atr_setup": {
    "prism_material": "ZnSe",      // Prism material name
    "prism_index": 2.4,            // Refractive index at measurement wavelength
    "incident_angle_deg": 60       // Incident angle in degrees
}
```

Common prism materials:
- ZnSe: n ≈ 2.4 (mid-IR)
- Ge: n ≈ 4.0 (mid-IR)
- Diamond: n ≈ 2.4 (mid-IR)

### Sample Parameters

For each sample, provide Hall effect measurements and file locations:

```json
"sample_A": {
    "label": "Descriptive name",
    "thickness_nm": 370,                      // Film thickness (nm)
    "carrier_concentration_cm3": 1.8e20,      // From Hall effect (cm⁻³)
    "mobility_cm2Vs": 10,                     // From Hall effect (cm²/V·s)
    "carrier_type": "p-type",                 // "p-type" or "n-type"
    "effective_mass_me": 0.8,                 // In units of electron mass
    "data_directory": "data/sample_A",        // Path to data files
    "data_pattern": "*.dpt",                  // Glob pattern to find files
    "fit_range_cm1": [300, 1000],            // Frequency range for fitting
    "plot_range_cm1": [200, 1200]            // Frequency range for plots
}
```

**Note:** The carrier concentration and mobility are used only for **initial guesses**. The final values are determined by fitting to the ATR data.

### Fitting Options

Control the optimization algorithm:

```json
"fitting_options": {
    "max_iterations": 500,              // Maximum optimization iterations
    "population_size": 15,              // Population size for DE algorithm
    "tolerance": 1e-6,                  // Convergence tolerance
    "random_seed": 42,                  // Random seed for reproducibility
    "omega_p_bounds_factor": [0.5, 2.0], // Search range for ω_p (× initial guess)
    "gamma_p_bounds_factor": [0.3, 3.0], // Search range for γ_p (× initial guess)
    "S_bounds": [0.5, 10.0],            // Oscillator strength bounds
    "gamma_ph_bounds": [10.0, 100.0]    // Phonon damping bounds (cm⁻¹)
}
```

### Output Options

```json
"output": {
    "directory": "output",          // Output directory
    "save_permittivity": true,      // Save permittivity text files
    "save_plots": true,             // Generate plots
    "plot_dpi": 150                 // Plot resolution
}
```

## Output Files

### Permittivity Data Files

Located in `output/permittivity_*.txt`:

```
# Extracted permittivity for Sample A
# ...fitted parameters...
# Columns: Frequency(cm^-1), Real(epsilon), Imag(epsilon)
400.000    10.234    -2.456
401.000    10.189    -2.478
...
```

These files can be directly imported into:
- Transfer Matrix Method (TMM) codes
- FDTD simulations (Lumerical, MEEP, etc.)
- Optical design software

### Visualization

`output/extraction_summary.png` contains:
1. **Left column:** Fit quality (experimental vs. fitted ATR reflectance)
2. **Middle column:** Real part of permittivity Re(ε)
3. **Right column:** Imaginary part of permittivity Im(ε)

## Method Details

### Fitting Algorithm

The code uses **differential evolution** (global optimization) to minimize the residual:

```
Residual = Σ[R_measured(ω) - R_calculated(ω)]²
```

For each trial parameter set:
1. Calculate ε(ω) using Drude-Lorentz model
2. Calculate ATR reflectance using Fresnel equations
3. Compare to experimental data

### Physical Interpretation

**Fitted parameters:**
- `ω_p`: Related to free carrier density: `ω_p ∝ √(n/m*)`
- `γ_p`: Related to carrier mobility: `γ_p ∝ 1/μ`
- `S`: Phonon oscillator strength (couples to electric field)
- `γ_ph`: Phonon damping (lattice quality, temperature effects)

**Fixed parameters:**
- `ε_∞`: High-frequency electronic contribution
- `ω_TO`: Lattice vibration frequency

### Assumptions

1. **Single Lorentz oscillator:** Appropriate for materials with one dominant phonon mode
2. **Thin film approximation:** Assumes simple Fresnel reflection (good for < 1 μm films)
3. **Isotropic material:** Doesn't account for anisotropy
4. **Uniform film:** No gradients or interfaces within the film

## Extending the Code

### Multiple Phonon Modes

To add additional Lorentz oscillators, modify the `drude_lorentz_permittivity()` function:

```python
# Add second phonon mode
lorentz2 = (S2 * omega_TO2**2) / (omega_TO2**2 - omega**2 - 1j * gamma_ph2 * omega)
eps = eps_inf + drude + lorentz1 + lorentz2
```

### Custom Reflectance Models

For more sophisticated geometries, modify `calculate_atr_reflectance()`:
- Multi-layer transfer matrix method
- Substrate effects
- Anisotropic materials

### Alternative Optimization

Replace `differential_evolution` with:
- `scipy.optimize.minimize` for local optimization
- Bayesian optimization for uncertainty quantification
- MCMC for full posterior distributions

## Troubleshooting

### "No files found for sample"

Check that:
- Data directory path is correct
- File pattern matches your file extension (.dpt, .txt, .csv)
- Files are in the expected location

### Poor fit quality

Try:
- Adjusting `fit_range_cm1` to focus on high-quality data regions
- Increasing `max_iterations` and `population_size`
- Widening parameter bounds if initial guesses are poor
- Checking for data artifacts (baseline offsets, noise)

### Unphysical results

- Verify Hall effect measurements (carrier concentration, mobility)
- Check `effective_mass_me` from literature
- Ensure `epsilon_inf` and `omega_TO` are appropriate for your material
- Consider whether single-oscillator model is appropriate

## Citation

If you use this code in your research, please cite:

Bartelsen, Emma. "ATR Permittivity Extraction Tool" (2025). 
GitHub repository: https://github.com/emmabartelsen123/Permittivity-Extraction-from-ATR-Measurements

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request

Areas for improvement:
- Support for multiple Lorentz oscillators
- Uncertainty quantification
- GUI for configuration
- Additional reflectance models (s-polarization, transmission)

## License

MIT License - see LICENSE file for details

**Emma Bartelsen**  
PhD Candidate, Interdisciplinary Materials Science  
Vanderbilt University  

Questions or issues? Please:
- Open a [GitHub issue](https://github.com/emmabartelsen123/Permittivity-Extraction-from-ATR-Measurements/issues)
- Email: emma.bartelsen@gmail.com

## References

1. Wooten, F. "Optical Properties of Solids" (Academic Press, 1972)
2. Palik, E.D. "Handbook of Optical Constants of Solids" (Academic Press, 1985)
3. Tompkins, H.G. & Irene, E.A. "Handbook of Ellipsometry" (William Andrew, 2005)

## Acknowledgments

Developed for characterizing doped semiconductor thin films for infrared photonics applications.
