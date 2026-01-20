#!/usr/bin/env python3
"""
Extract Complex Permittivity ε(ω) from ATR Reflectance Spectra
Uses Drude-Lorentz Model Fitting

For doped semiconductors:
ε(ω) = ε_∞ - (ω_p² / (ω² + iγ_p·ω)) + Σ[(S_j·ω_TO,j²) / (ω_TO,j² - ω² - iγ_j·ω)]

Where:
- Drude term: Free carrier contribution (plasmon)
- Lorentz terms: Phonon contributions

Author: [Your Name]
License: MIT
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution
from scipy.constants import e, epsilon_0, m_e
import os
import glob
import json

# ============================================================================
# USER CONFIGURATION
# ============================================================================

def load_config(config_file='config.json'):
    """
    Load configuration from JSON file
    
    If config file doesn't exist, creates a template and exits.
    """
    if not os.path.exists(config_file):
        print(f"\n⚠️  Configuration file '{config_file}' not found!")
        print("Creating template config file...\n")
        create_config_template(config_file)
        print(f"✅ Template created: {config_file}")
        print("\nPlease edit this file with your sample parameters and run again.\n")
        exit()
    
    with open(config_file, 'r') as f:
        config = json.load(f)
    
    return config

def create_config_template(config_file='config.json'):
    """Create a template configuration file"""
    template = {
        "material_constants": {
            "name": "ScN",
            "epsilon_inf": 12.8,
            "omega_TO": 359.0,
            "omega_LO": 677.0,
            "description": "High-frequency dielectric constant and phonon frequencies (cm^-1)"
        },
        "atr_setup": {
            "prism_material": "ZnSe",
            "prism_index": 2.4,
            "incident_angle_deg": 60,
            "description": "ATR experimental configuration"
        },
        "samples": {
            "sample_A": {
                "label": "Sample A - High-doped p-type",
                "thickness_nm": 370,
                "carrier_concentration_cm3": 1.8e20,
                "mobility_cm2Vs": 10,
                "carrier_type": "p-type",
                "effective_mass_me": 0.8,
                "data_directory": "data/sample_A",
                "data_pattern": "*.dpt",
                "fit_range_cm1": [300, 1000],
                "plot_range_cm1": [200, 1200]
            },
            "sample_B": {
                "label": "Sample B - Example n-type",
                "thickness_nm": 135,
                "carrier_concentration_cm3": 6.8e20,
                "mobility_cm2Vs": 73.9,
                "carrier_type": "n-type",
                "effective_mass_me": 0.39,
                "data_directory": "data/sample_B",
                "data_pattern": "*.dpt",
                "fit_range_cm1": [2000, 5000],
                "plot_range_cm1": [500, 5000]
            }
        },
        "fitting_options": {
            "max_iterations": 500,
            "population_size": 15,
            "tolerance": 1e-6,
            "random_seed": 42,
            "omega_p_bounds_factor": [0.5, 2.0],
            "gamma_p_bounds_factor": [0.3, 3.0],
            "S_bounds": [0.5, 10.0],
            "gamma_ph_bounds": [10.0, 100.0],
            "description": "Optimization parameters for differential evolution"
        },
        "output": {
            "directory": "output",
            "save_permittivity": true,
            "save_plots": true,
            "plot_dpi": 150
        }
    }
    
    with open(config_file, 'w') as f:
        json.dump(template, f, indent=4)

# ============================================================================
# DRUDE-LORENTZ MODEL
# ============================================================================

def cm_to_angular_freq(wavenumber_cm):
    """Convert wavenumber (cm^-1) to angular frequency (rad/s)"""
    c = 2.998e10  # cm/s
    return 2 * np.pi * c * wavenumber_cm

def drude_lorentz_permittivity(omega_cm, params):
    """
    Calculate complex permittivity using Drude-Lorentz model
    
    Parameters:
    -----------
    omega_cm : array, wavenumber in cm^-1
    params : dict with keys
        'eps_inf': high-frequency dielectric constant
        'omega_p': plasma frequency (cm^-1)
        'gamma_p': plasma damping (cm^-1)
        'S': oscillator strength for phonon
        'omega_TO': TO phonon frequency (cm^-1)
        'gamma_ph': phonon damping (cm^-1)
    
    Returns:
    --------
    eps : complex array, complex permittivity
    """
    omega = omega_cm
    eps_inf = params['eps_inf']
    omega_p = params['omega_p']
    gamma_p = params['gamma_p']
    S = params['S']
    omega_TO = params['omega_TO']
    gamma_ph = params['gamma_ph']
    
    # Drude term (free carriers)
    drude = -omega_p**2 / (omega**2 + 1j * gamma_p * omega)
    
    # Lorentz term (phonon)
    lorentz = (S * omega_TO**2) / (omega_TO**2 - omega**2 - 1j * gamma_ph * omega)
    
    eps = eps_inf + drude + lorentz
    
    return eps

def calculate_plasma_frequency(n_cm3, m_eff):
    """
    Calculate plasma frequency from carrier concentration
    
    omega_p = sqrt(n * e^2 / (epsilon_0 * m*))
    
    Parameters:
    -----------
    n_cm3 : float, carrier concentration in cm^-3
    m_eff : float, effective mass in units of m_e
    
    Returns:
    --------
    omega_p_cm : float, plasma frequency in cm^-1
    """
    n = n_cm3 * 1e6  # convert cm^-3 to m^-3
    m_eff_kg = m_eff * m_e
    omega_p_rad = np.sqrt(n * e**2 / (epsilon_0 * m_eff_kg))  # rad/s
    c = 2.998e10  # cm/s
    omega_p_cm = omega_p_rad / (2 * np.pi * c)
    return omega_p_cm

def calculate_damping_from_mobility(mu_cm2Vs, m_eff):
    """
    Calculate damping rate from mobility
    
    gamma = e / (mu * m*)
    
    Parameters:
    -----------
    mu_cm2Vs : float, mobility in cm²/V-s
    m_eff : float, effective mass in units of m_e
    
    Returns:
    --------
    gamma_cm : float, damping rate in cm^-1
    """
    mu = mu_cm2Vs * 1e-4  # convert cm²/V-s to m²/V-s
    m_eff_kg = m_eff * m_e
    gamma_rad = e / (mu * m_eff_kg)  # rad/s
    c = 2.998e10  # cm/s
    gamma_cm = gamma_rad / (2 * np.pi * c)
    return gamma_cm

# ============================================================================
# ATR REFLECTANCE CALCULATION
# ============================================================================

def calculate_atr_reflectance(omega_cm, eps_film, theta_deg=60, n_prism=2.4):
    """
    Calculate ATR reflectance for prism/film/substrate geometry
    
    Simplified model: Fresnel reflection at prism/film interface
    
    Parameters:
    -----------
    omega_cm : wavenumber (cm^-1)
    eps_film : complex permittivity of film
    theta_deg : incident angle in prism (degrees)
    n_prism : refractive index of prism
    
    Returns:
    --------
    R : reflectance (p-polarization)
    """
    theta = np.radians(theta_deg)
    
    # Wave vector components
    k0 = omega_cm  # in units of cm^-1
    kx = n_prism * k0 * np.sin(theta)
    
    # Perpendicular components
    kz_prism = np.sqrt((n_prism * k0)**2 - kx**2 + 0j)
    kz_film = np.sqrt(eps_film * k0**2 - kx**2 + 0j)
    
    # Fresnel coefficient for p-polarization (TM mode)
    r_p = (eps_film * kz_prism - n_prism**2 * kz_film) / \
          (eps_film * kz_prism + n_prism**2 * kz_film)
    
    R = np.abs(r_p)**2
    
    return R

# ============================================================================
# FITTING ROUTINE
# ============================================================================

def objective_function(params_array, omega_data, R_data, material_constants, 
                      atr_config):
    """
    Objective function for fitting: sum of squared residuals
    
    Parameters:
    -----------
    params_array : array, [omega_p, gamma_p, S, gamma_ph]
    omega_data : array, experimental wavenumbers
    R_data : array, experimental reflectance
    material_constants : dict, material parameters
    atr_config : dict, ATR setup parameters
    """
    # Unpack parameters
    omega_p, gamma_p, S, gamma_ph = params_array
    
    params = {
        'eps_inf': material_constants['epsilon_inf'],
        'omega_p': omega_p,
        'gamma_p': gamma_p,
        'S': S,
        'omega_TO': material_constants['omega_TO'],
        'gamma_ph': gamma_ph
    }
    
    # Calculate permittivity
    eps = drude_lorentz_permittivity(omega_data, params)
    
    # Calculate reflectance
    R_calc = calculate_atr_reflectance(
        omega_data, eps, 
        theta_deg=atr_config['incident_angle_deg'],
        n_prism=atr_config['prism_index']
    )
    
    # Residual
    residual = np.sum((R_calc - R_data)**2)
    
    return residual

def fit_permittivity_to_atr(omega_data, R_data, sample_info, material_constants,
                            atr_config, fit_options):
    """
    Fit Drude-Lorentz parameters to ATR reflectance data
    
    Parameters:
    -----------
    omega_data : array, experimental wavenumbers
    R_data : array, experimental reflectance
    sample_info : dict, sample parameters
    material_constants : dict, material parameters
    atr_config : dict, ATR setup
    fit_options : dict, optimization parameters
    
    Returns:
    --------
    fitted_params : dict with fitted parameters
    """
    print(f"\n{'='*70}")
    print(f"FITTING: {sample_info['label']}")
    print(f"{'='*70}")
    
    # Initial guesses from sample parameters
    omega_p_init = calculate_plasma_frequency(
        sample_info['carrier_concentration_cm3'], 
        sample_info['effective_mass_me']
    )
    gamma_p_init = calculate_damping_from_mobility(
        sample_info['mobility_cm2Vs'],
        sample_info['effective_mass_me']
    )
    
    print(f"Initial guesses from Hall measurements:")
    print(f"  Plasma frequency: {omega_p_init:.1f} cm⁻¹")
    print(f"  Plasma damping: {gamma_p_init:.1f} cm⁻¹")
    
    # Parameter bounds: [omega_p, gamma_p, S, gamma_ph]
    bounds = [
        (omega_p_init * fit_options['omega_p_bounds_factor'][0], 
         omega_p_init * fit_options['omega_p_bounds_factor'][1]),
        (gamma_p_init * fit_options['gamma_p_bounds_factor'][0], 
         gamma_p_init * fit_options['gamma_p_bounds_factor'][1]),
        tuple(fit_options['S_bounds']),
        tuple(fit_options['gamma_ph_bounds']),
    ]
    
    print(f"\nRunning global optimization (differential_evolution)...")
    
    # Global optimization
    result = differential_evolution(
        objective_function,
        bounds,
        args=(omega_data, R_data, material_constants, atr_config),
        maxiter=fit_options['max_iterations'],
        popsize=fit_options['population_size'],
        tol=fit_options['tolerance'],
        seed=fit_options['random_seed'],
        workers=1,
        updating='deferred',
        disp=True
    )
    
    print(f"\nOptimization complete!")
    print(f"  Final residual: {result.fun:.6e}")
    
    # Extract fitted parameters
    omega_p_fit, gamma_p_fit, S_fit, gamma_ph_fit = result.x
    
    fitted_params = {
        'eps_inf': material_constants['epsilon_inf'],
        'omega_p': omega_p_fit,
        'gamma_p': gamma_p_fit,
        'S': S_fit,
        'omega_TO': material_constants['omega_TO'],
        'gamma_ph': gamma_ph_fit
    }
    
    print(f"\n{'='*70}")
    print(f"FITTED PARAMETERS:")
    print(f"{'='*70}")
    print(f"  ε_∞ = {material_constants['epsilon_inf']:.2f} (fixed)")
    print(f"  ω_p = {omega_p_fit:.1f} cm⁻¹ (plasma frequency)")
    print(f"  γ_p = {gamma_p_fit:.1f} cm⁻¹ (plasma damping)")
    print(f"  S = {S_fit:.2f} (oscillator strength)")
    print(f"  ω_TO = {material_constants['omega_TO']:.1f} cm⁻¹ (fixed)")
    print(f"  γ_ph = {gamma_ph_fit:.1f} cm⁻¹ (phonon damping)")
    print(f"{'='*70}\n")
    
    return fitted_params

# ============================================================================
# DATA I/O
# ============================================================================

def load_atr_data(filepath, delimiter=','):
    """
    Load ATR data file
    
    Expected format: wavenumber, reflectance (comma or tab separated)
    """
    try:
        data = np.loadtxt(filepath, delimiter=delimiter)
    except:
        # Try tab delimiter if comma fails
        data = np.loadtxt(filepath, delimiter='\t')
    
    wavenumber = data[:, 0]
    reflectance = data[:, 1]
    return wavenumber, reflectance

def save_permittivity(omega, eps, fitted_params, sample_info, output_dir, sample_id):
    """Save extracted permittivity to file"""
    filename = os.path.join(output_dir, f'permittivity_{sample_id}.txt')
    
    header = (
        f"Extracted permittivity for {sample_info['label']}\n"
        f"Material: {sample_info.get('material', 'N/A')}\n"
        f"Thickness: {sample_info['thickness_nm']:.1f} nm\n"
        f"Carrier concentration: {sample_info['carrier_concentration_cm3']:.2e} cm^-3\n"
        f"Mobility: {sample_info['mobility_cm2Vs']:.2f} cm²/V-s\n"
        f"Carrier type: {sample_info['carrier_type']}\n"
        f"\nFitted Drude-Lorentz parameters:\n"
        f"  omega_p = {fitted_params['omega_p']:.2f} cm^-1 (plasma frequency)\n"
        f"  gamma_p = {fitted_params['gamma_p']:.2f} cm^-1 (plasma damping)\n"
        f"  S = {fitted_params['S']:.2f} (oscillator strength)\n"
        f"  gamma_ph = {fitted_params['gamma_ph']:.2f} cm^-1 (phonon damping)\n"
        f"  eps_inf = {fitted_params['eps_inf']:.2f} (fixed)\n"
        f"  omega_TO = {fitted_params['omega_TO']:.2f} cm^-1 (fixed)\n"
        f"\nColumns: Frequency(cm^-1), Real(epsilon), Imag(epsilon)"
    )
    
    np.savetxt(filename, 
               np.column_stack([omega, eps.real, eps.imag]),
               header=header,
               fmt='%.6f',
               delimiter='\t')
    
    return filename

# ============================================================================
# PLOTTING
# ============================================================================

def plot_results(all_results, material_constants, output_dir, dpi=150):
    """Generate comprehensive plots of fitting results"""
    
    n_samples = len(all_results)
    fig = plt.figure(figsize=(20, 4*n_samples))
    
    for idx, (sample_id, results) in enumerate(all_results.items()):
        
        # Plot 1: ATR fit quality
        ax1 = plt.subplot(n_samples, 3, idx*3 + 1)
        ax1.plot(results['omega_data'], results['R_data'], 'o', 
                markersize=3, alpha=0.5, label='Experimental', color='steelblue')
        ax1.plot(results['omega_data'], results['R_fitted'], '-', 
                linewidth=2, label='Fitted', color='crimson')
        ax1.set_xlabel('Wavenumber (cm⁻¹)', fontsize=11)
        ax1.set_ylabel('Reflectance', fontsize=11)
        ax1.set_title(f'{sample_id} - ATR Fit Quality', fontweight='bold')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Plot 2: Real part of permittivity
        ax2 = plt.subplot(n_samples, 3, idx*3 + 2)
        ax2.plot(results['omega'], results['epsilon'].real, '-', 
                linewidth=2, color='navy')
        ax2.axhline(0, color='k', linestyle='--', alpha=0.3)
        if 'omega_TO' in material_constants:
            ax2.axvline(material_constants['omega_TO'], color='r', 
                       linestyle='--', alpha=0.5, label='ω_TO')
        if 'omega_LO' in material_constants:
            ax2.axvline(material_constants['omega_LO'], color='g', 
                       linestyle='--', alpha=0.5, label='ω_LO')
        ax2.set_xlabel('Wavenumber (cm⁻¹)', fontsize=11)
        ax2.set_ylabel('Re(ε)', fontsize=11)
        ax2.set_title(f'{sample_id} - Real Permittivity', fontweight='bold')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # Plot 3: Imaginary part of permittivity
        ax3 = plt.subplot(n_samples, 3, idx*3 + 3)
        ax3.plot(results['omega'], results['epsilon'].imag, '-', 
                linewidth=2, color='darkred')
        if 'omega_TO' in material_constants:
            ax3.axvline(material_constants['omega_TO'], color='r', 
                       linestyle='--', alpha=0.5, label='ω_TO')
        if 'omega_LO' in material_constants:
            ax3.axvline(material_constants['omega_LO'], color='g', 
                       linestyle='--', alpha=0.5, label='ω_LO')
        ax3.set_xlabel('Wavenumber (cm⁻¹)', fontsize=11)
        ax3.set_ylabel('Im(ε)', fontsize=11)
        ax3.set_title(f'{sample_id} - Imaginary Permittivity', fontweight='bold')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plot_filename = os.path.join(output_dir, 'extraction_summary.png')
    plt.savefig(plot_filename, dpi=dpi, bbox_inches='tight')
    plt.close()
    
    return plot_filename

# ============================================================================
# MAIN EXTRACTION ROUTINE
# ============================================================================

def main(config_file='config.json'):
    """Main extraction routine"""
    
    print("\n" + "="*80)
    print("PERMITTIVITY EXTRACTION FROM ATR REFLECTANCE")
    print("="*80)
    
    # Load configuration
    config = load_config(config_file)
    
    material_constants = config['material_constants']
    atr_config = config['atr_setup']
    fit_options = config['fitting_options']
    output_config = config['output']
    
    # Create output directory
    output_dir = output_config['directory']
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    all_results = {}
    
    # Process each sample
    for sample_id, sample_info in config['samples'].items():
        
        print(f"\n{'='*80}")
        print(f"PROCESSING: {sample_id}")
        print(f"{'='*80}")
        
        # Find data files
        data_dir = sample_info['data_directory']
        pattern = sample_info['data_pattern']
        file_pattern = os.path.join(data_dir, pattern)
        files = glob.glob(file_pattern)
        
        if not files:
            print(f"\n⚠️  WARNING: No files found for {sample_id}")
            print(f"   Pattern: {file_pattern}")
            print(f"   Skipping...\n")
            continue
        
        filepath = files[0]
        print(f"\nData file: {os.path.basename(filepath)}")
        
        # Load ATR data
        omega_data, R_data = load_atr_data(filepath)
        
        # Apply frequency range for fitting
        fit_range = sample_info['fit_range_cm1']
        mask = (omega_data >= fit_range[0]) & (omega_data <= fit_range[1])
        omega_fit = omega_data[mask]
        R_fit = R_data[mask]
        
        print(f"Data points for fitting: {len(omega_fit)}")
        print(f"Frequency range: {omega_fit.min():.0f} - {omega_fit.max():.0f} cm⁻¹")
        
        # Fit permittivity
        fitted_params = fit_permittivity_to_atr(
            omega_fit, R_fit, sample_info, 
            material_constants, atr_config, fit_options
        )
        
        # Generate full permittivity curve
        plot_range = sample_info['plot_range_cm1']
        omega_full = np.linspace(plot_range[0], plot_range[1], 2000)
        eps_full = drude_lorentz_permittivity(omega_full, fitted_params)
        
        # Calculate fitted reflectance
        R_fitted = calculate_atr_reflectance(
            omega_fit, 
            drude_lorentz_permittivity(omega_fit, fitted_params),
            theta_deg=atr_config['incident_angle_deg'],
            n_prism=atr_config['prism_index']
        )
        
        # Store results
        all_results[sample_id] = {
            'params': fitted_params,
            'omega': omega_full,
            'epsilon': eps_full,
            'omega_data': omega_fit,
            'R_data': R_fit,
            'R_fitted': R_fitted,
            'info': sample_info
        }
        
        # Save permittivity data
        if output_config['save_permittivity']:
            eps_filename = save_permittivity(
                omega_full, eps_full, fitted_params, 
                sample_info, output_dir, sample_id
            )
            print(f"✅ Saved: {eps_filename}")
    
    # Generate plots
    if output_config['save_plots'] and all_results:
        print("\n" + "="*80)
        print("GENERATING PLOTS")
        print("="*80)
        
        plot_filename = plot_results(
            all_results, material_constants, output_dir, 
            dpi=output_config['plot_dpi']
        )
        print(f"✅ Saved plot: {plot_filename}")
    
    # Print summary
    print("\n" + "="*80)
    print("EXTRACTION SUMMARY")
    print("="*80)
    
    for sample_id, results in all_results.items():
        params = results['params']
        print(f"\n{sample_id}: {results['info']['label']}")
        print(f"  ω_p = {params['omega_p']:.1f} cm⁻¹ (plasma)")
        print(f"  γ_p = {params['gamma_p']:.1f} cm⁻¹ (plasma damping)")
        print(f"  S = {params['S']:.2f} (oscillator strength)")
        print(f"  γ_ph = {params['gamma_ph']:.1f} cm⁻¹ (phonon damping)")
    
    print("\n" + "="*80)
    print("✅ PERMITTIVITY EXTRACTION COMPLETE!")
    print("="*80)
    print(f"\nOutput files saved to: {output_dir}/")
    print(f"\nExtracted permittivity files can be used in:")
    print(f"  • Transfer Matrix Method (TMM) simulations")
    print(f"  • FDTD electromagnetic simulations")
    print(f"  • Optical modeling and design")
    print("="*80 + "\n")

if __name__ == "__main__":
    import sys
    
    # Allow custom config file via command line
    config_file = sys.argv[1] if len(sys.argv) > 1 else 'config.json'
    main(config_file)
