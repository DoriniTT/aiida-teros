"""
Configuration schema for AiiDA-TEROS.

This module defines the JSON schema for validating configuration files.
"""

CONFIG_SCHEMA = {
    "$schema": "http://json-schema.org/draft-07/schema#",
    "type": "object",
    "required": [
        "bulk_structure_path",
        "bulk_metal_structure_path",
        "potential_family",
        "code_label",
        "thermodynamic_parameters",
        "incar_parameters_bulk",
        "incar_parameters_bulk_metal",
        "incar_parameters_o2",
        "incar_parameters_slab",
        "workflow_settings",
        "potential_mapping",
        "parser_settings",
        "computer_options",
        "slab_parameters"
    ],
    "properties": {
        "bulk_structure_path": {
            "type": "string",
            "description": "Path to the bulk structure VASP file"
        },
        "bulk_metal_structure_path": {
            "type": "string",
            "description": "Path to the bulk metal structure VASP file"
        },
        "potential_family": {
            "type": "string",
            "description": "Potential family for VASP calculations"
        },
        "code_label": {
            "type": "string",
            "description": "Label of the VASP code configured in AiiDA"
        },
        "thermodynamic_parameters": {
            "type": "object",
            "required": ["hf_bulk"],
            "properties": {
                "hf_bulk": {
                    "type": "number",
                    "description": "Heat of formation for bulk in eV per formula unit"
                }
            }
        },
        "total_energies": {
            "type": "object",
            "properties": {
                "total_energy_first_element": {
                    "type": "number",
                    "description": "Energy of the first element (e.g., Ag) in eV - now calculated automatically"
                },
                "total_energy_o2": {
                    "type": "number",
                    "description": "Energy for O2 molecule in eV - now calculated automatically"
                }
            }
        },
        "terminations": {
            "type": "object",
            "description": "Optional specific terminations to study",
            "additionalProperties": {
                "type": "string",
                "description": "Path to the termination VASP file"
            }
        },
        "incar_parameters_bulk": {
            "type": "object",
            "required": ["incar"],
            "properties": {
                "incar": {
                    "type": "object",
                    "description": "INCAR parameters for bulk relaxation"
                }
            }
        },
        "incar_parameters_bulk_metal": {
            "type": "object",
            "required": ["incar"],
            "properties": {
                "incar": {
                    "type": "object",
                    "description": "INCAR parameters for bulk metal relaxation"
                }
            }
        },
        "incar_parameters_o2": {
            "type": "object",
            "required": ["incar"],
            "properties": {
                "incar": {
                    "type": "object",
                    "description": "INCAR parameters for O2 molecule relaxation"
                }
            }
        },
        "incar_parameters_slab": {
            "type": "object",
            "required": ["incar"],
            "properties": {
                "incar": {
                    "type": "object",
                    "description": "INCAR parameters for slab relaxation"
                }
            }
        },
        "workflow_settings": {
            "type": "object",
            "required": ["kpoints_precision"],
            "properties": {
                "kpoints_precision": {
                    "type": "number",
                    "description": "K-points mesh density precision"
                }
            }
        },
        "potential_mapping": {
            "type": "object",
            "description": "Mapping of potentials for each element"
        },
        "parser_settings": {
            "type": "object",
            "description": "Settings for the VASP parser"
        },
        "computer_options": {
            "type": "object",
            "required": ["resources"],
            "properties": {
                "resources": {
                    "type": "object",
                    "required": ["num_machines", "num_cores_per_machine"],
                    "properties": {
                        "num_machines": {
                            "type": "integer",
                            "description": "Number of machines for calculation"
                        },
                        "num_cores_per_machine": {
                            "type": "integer",
                            "description": "Number of cores per machine"
                        }
                    }
                },
                "queue_name": {
                    "type": "string",
                    "description": "Queue name for job submission"
                }
            }
        },
        "slab_parameters": {
            "type": "object",
            "required": ["miller_indices", "min_slab_thickness"],
            "properties": {
                "miller_indices": {
                    "type": "array",
                    "items": {"type": "integer"},
                    "minItems": 3,
                    "maxItems": 3,
                    "description": "Miller indices for slab generation"
                },
                "min_slab_thickness": {
                    "type": "number",
                    "description": "Minimum slab thickness in Angstroms"
                }
            }
        }
    }
}