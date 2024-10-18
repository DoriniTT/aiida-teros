# MasterThermoWorkChain: Automated Surface Thermodynamics with AiiDA-VASP

Welcome to the **MasterThermoWorkChain** repository! This project provides a comprehensive workflow for automating surface thermodynamics calculations using AiiDA and VASP. It allows users to:

- Generate symmetric surface terminations from a given bulk structure and orientation.
- Perform relaxation calculations on all generated slabs using VASP.
- Compute surface Gibbs free energies and construct surface phase diagrams using ab initio atomistic thermodynamics.
- Identify and output the most stable surface terminations based on thermodynamic analysis.

This README will guide you through the prerequisites, installation, and usage of the MasterThermoWorkChain.

---

## Table of Contents

- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
  - [1. Setting Up the Environment](#1-setting-up-the-environment)
  - [2. Preparing Input Files](#2-preparing-input-files)
  - [3. Running the WorkChain](#3-running-the-workchain)
  - [4. Monitoring and Retrieving Results](#4-monitoring-and-retrieving-results)
- [Output Explanation](#output-explanation)
- [Examples](#examples)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

---

## Prerequisites

Before using the MasterThermoWorkChain, ensure you have the following software and packages installed:

- **Python 3.7+**
- **AiiDA** (v1.0 or higher)
- **VASP** (and a valid license)
- **AiiDA-VASP** plugin
- **ASE** (Atomic Simulation Environment)
- **pymatgen** (Python Materials Genomics)
- **NumPy**, **Matplotlib**, **Seaborn**
- **pint** (for unit handling)
- **Additional Python packages**: `subprocess`, `shutil`, `os`, `time`

### Note on VASP

VASP is a proprietary software and requires a valid license to use. Ensure that you have access to VASP and that it is correctly set up with AiiDA on your computing resources.

---

## Installation

1. **Clone the Repository**

   ```bash
   git clone https://github.com/your_username/MasterThermoWorkChain.git
   cd MasterThermoWorkChain
