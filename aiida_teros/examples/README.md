# AiiDA-TEROS Examples

This directory contains examples demonstrating how to use AiiDA-TEROS for surface thermodynamics calculations of different oxide materials.

## Available Examples

1. **[Ag2O](./ag2o/)** - Binary oxide example
   - Demonstrates calculations for a simple Ag2O binary oxide
   - Shows how to set up and run binary oxide calculations
   - Includes example structure files and configuration

2. **[Ag3PO4](./ag3po4/)** - Ternary oxide example
   - Demonstrates calculations for a more complex Ag3PO4 ternary oxide
   - Shows how to set up and run ternary oxide calculations
   - Includes example structure files and configuration

## How to Use These Examples

Each example directory contains:
- A comprehensive README.md with step-by-step instructions
- Structure files in VASP POSCAR format
- Configuration file (config.yaml)
- Python script to run the example (run_example.py)

To run an example:

1. Navigate to the example directory:
   ```bash
   cd ag2o/  # or cd ag3po4/
   ```

2. Review and modify the configuration file as needed:
   ```bash
   nano config.yaml
   ```

3. Run the example:
   ```bash
   python run_example.py
   ```

## Common Configuration Parameters

All examples share a similar configuration structure:

```yaml
# Paths and File Names
bulk_structure_path: "oxide.vasp"            # Path to bulk oxide structure
bulk_metal_structure_path: "metal.vasp"      # Path to bulk metal structure
potential_family: "PBE"                      # Potential family
code_label: "vasp@computer"                  # VASP code label in AiiDA

# Thermodynamic Parameters
thermodynamic_parameters:
  hf_bulk: -X.XX                             # Heat of formation for bulk in eV

# Slab Generation Parameters
slab_parameters:
  miller_indices: [h, k, l]                  # Miller indices
  min_slab_thickness: 10.0                   # Minimum slab thickness in Å
  vacuum: 15.0                               # Vacuum thickness in Å
```

## Extending the Examples

You can modify these examples to study your own materials:

1. Replace the structure files with your own materials
2. Adjust the configuration parameters as needed
3. Run the workflow following the same process

## Additional Resources

For more information on AiiDA-TEROS, refer to:
- The main AiiDA-TEROS documentation
- The repository README.md
- The API documentation