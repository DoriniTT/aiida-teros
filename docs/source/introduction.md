# Introduction

**AiiDA-TEROS** stands for **Thermodynamics of Oxide Surfaces**, and this project is dedicated to automating the study of surface thermodynamics for oxide materials using the powerful combination of [AiiDA](https://www.aiida.net/) and [VASP](https://www.vasp.at/).

## Why Use AiiDA-TEROS?

Understanding the thermodynamic properties of oxide surfaces is crucial in fields such as catalysis, corrosion science, and materials engineering. Oxide surfaces often exhibit complex behaviors due to their interactions with the environment, especially in the presence of oxygen and under varying temperature and pressure conditions. Traditional computational studies of these systems can be time-consuming and labor-intensive, requiring significant manual setup and data handling.

AiiDA-TEROS addresses these challenges by providing:

- **Efficiency**: Automate repetitive and complex tasks involved in surface thermodynamics calculations
- **Scalability**: Easily extend the workflow to study various oxide materials
- **Reproducibility**: Ensure consistent computational procedures
- **Customization**: Tailor calculation parameters to fit specific research needs

## Key Features

- **Automated Surface Generation**: Generate symmetric surface terminations from any bulk structure and crystallographic orientation
- **Efficient Relaxation Calculations**: Perform relaxation calculations on all generated slabs using VASP
- **Surface Thermodynamics Analysis**: Compute surface Gibbs free energies and construct surface phase diagrams
- **Stable Surface Identification**: Automatically identify the most thermodynamically stable surface terminations

![Workflow diagram](../images/fluxogram.pdf)

## Architecture

AiiDA-TEROS follows a modular architecture to improve maintainability and extensibility. The main components include:

- **Core Components**: Base classes and abstractions
- **Workflows**: AiiDA workflow implementations for surface thermodynamics
- **Utilities**: Helper functions for structure manipulation, thermodynamics, etc.
- **Schemas**: Data validation using JSON schemas
- **Plotting**: Visualization tools for surface energies and phase diagrams

By automating these complex tasks, AiiDA-TEROS accelerates the research process, enabling you to explore new materials and surface phenomena more efficiently.