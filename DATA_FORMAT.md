# Example ATR Data File Format

## Format

ATR data files should contain two columns:
1. Wavenumber (cm⁻¹)
2. Reflectance (0-1 or 0-100%)

Supported delimiters: comma (,) or tab

## Example .dpt file content:

```
400.0, 0.234
401.0, 0.236
402.0, 0.238
403.0, 0.241
404.0, 0.245
...
```

Or with tab delimiter:

```
400.0	0.234
401.0	0.236
402.0	0.238
403.0	0.241
404.0	0.245
...
```

## File extensions

Common extensions that work:
- `.dpt` (typical FTIR output)
- `.txt`
- `.csv`
- `.dat`

The code will automatically detect the delimiter.

## Generating example data

If you want to test the code without real ATR data, you can generate synthetic data:

```python
import numpy as np

# Synthetic ATR spectrum
omega = np.linspace(300, 1000, 500)
R = 0.3 + 0.1 * np.exp(-((omega - 600)/50)**2)  # Gaussian dip
R += 0.02 * np.random.randn(len(omega))  # Add noise

# Save
np.savetxt('example_data.dpt', np.column_stack([omega, R]), 
           delimiter=',', fmt='%.2f')
```

This creates a reflectance spectrum with a phonon-like resonance feature near 600 cm⁻¹.
