# Molecular Docking and Dynamics for E.COLLIPHA project

This repository corresponds to one of the three computational biology analyses of the E.COLLIPHA project in the GOGEC Competition of 2026.

> **Objective**: Identify the potential mechanical properties of PHA as a polymer in medical sutures for wound healing.

We worked using VS Code, along with the WSL extension (Ubuntu 22.04) as our development environment (Windows users only).

1. WSL Installation: https://learn.microsoft.com/en-us/windows/wsl/install
2. VS Code Installation: https://code.visualstudio.com/docs/setup/windows#_install-vs-code-on-windows
3. WSL extension on VS Code: https://code.visualstudio.com/docs/remote/wsl

- E.COLLIPHA Members involved (No Specific Order): Hannia Jeudi Martínez Vázquez, Victoria Alejandra Saucedo Farias

- Advisor: LBG Mauricio Reyes-Elizondo 

- PI: Dr J Claudio Moreno-Rocha

## Activity 0: Setting-Up Python Environment

We set a series of tutorials to create Conda-based Python environment in ```docking_tutorials/```. With these tutorials, we were able to work with a stable version of ```Autodock Vina``` and collaborate remotely using Git.

Run notebooks sequentially (Activity 0 → Activity 1 → Activity 2 → Activity 3):

```bash
   jupyter notebook 00_minimized_pha_oligomer.ipynb
   jupyter notebook 001_minimized_pga_validation.ipynb
   jupyter notebook 01_docking_human_proteins.ipynb
   jupyter notebook 002_docking_pga_validation.ipynb
   jupyter notebook 02_docking_results_analysis.ipynb
   jupyter notebook 03_mechanical_properties.ipynb
```

## Activity 1: Energy-minimize PHB oligomers

**Notebook**: ```00_minimized_pha_oligomer.ipynb``` and ```001_minimized_pga_validation.ipynb```

### Objective
> Generate energy-minimized 3D structures of PHB (polyhydroxybutyrate) oligomers ranging from 1 to 12 units for use as ligands in molecular docking. Energy-minimized PGA oligomers (1-12 units) as validation.

### Methodology

The PHB monomer SMILES was retrieved from PubChem: `C[C@H](CC(=O)O)O`. The polymerization reaction was defined using SMARTS notation to model esterification between carboxylic acid and alcohol functional groups:

```
Reaction SMARTS: [C:1](=[O:2])[OH].[C:3][OH:4]>>[C:1](=[O:2])[O:4][C:3]
```

PHB oligomers with 1 to 12 units were iteratively synthesized using ```RDKit``` by repeatedly applying the polymerization reaction. Each oligomer was sanitized and validated after each reaction step.

Three-dimensional coordinates were generated for each oligomer through the following pipeline:

- Addition of all hydrogen atoms using `Chem.AddHs()`
- 3D coordinate embedding using ETKDG algorithm (`AllChem.EmbedMolecule()`)
- Preliminary geometry optimization using MMFF94 force field (`AllChem.MMFFOptimizeMolecule()`)

Energy minimization was performed using ```OpenBabel``` with the MMFF94 force field for 500 steps (`localopt()` with `forcefield='mmff94', steps=500`). This step refined the 3D geometries to represent stable conformations.

All energy-minimized oligomer structures were exported to PDB format for downstream docking studies, stored in ```output/ligand_files/```.

## Activity 2: Docking of PHA oligomers with skin and immune proteins

**Notebook**: ```01_docking_human_proteins.ipynb``` and ```002_docking_pga_validation.ipynb```

### Objective
> Perform blind molecular docking of PHA oligomers (1-12 units) against eight human proteins involved in immune response and skin structure to assess binding interactions. The results were compared using PGA oligomers (1-12 units) as validation.

### Methodology

Energy-minimized PDB files from Activity 1 were converted to PDBQT format using ```OpenBabel``` with the following parameters:

- pH: ```7.4``` (physiologically relevant)
- Charge method: ```Gasteiger``` (semi-empirical partial charges)
- Hydrogen treatment: ```Polar hydrogens retained for accurate electrostatics```

Eight human protein structures were selected for docking analysis:

| PDB ID | Protein Name | Function |
|--------|--------------|----------|
| 6EC0 | Keratin 1 | Structural protein (skin) |
| 7CWK | Collagen type I | Extracellular matrix |
| 4BSO | R-SPONDIN-1 | Growth signaling |
| 1RG8 | Heparin-binding growth factor 1 | Cell growth/differentiation |
| 1CSG | Granulocyte-macrophage colony-stimulating factor | Immune response |
| 1TGJ | Transforming growth factor-beta 3 | Wound healing/immune regulation |
| 1VPF | Vascular endothelial growth factor | Angiogenesis |
| 2Z80 | Toll-like receptor 2 | Pattern recognition (immune) |

PDB structures were downloaded from the PDB database using the `requests` library. Downloaded PDB files were pre-processed to:

- Remove non-aminoacidic residues (ligands, water molecules, ions)
- Retain protein chains while removing duplicate structures
- Remove polar hydrogens
- Standardize pH to 7.4
- Apply Gasteiger charge assignment

The processed protein structures were converted to PDBQT format using ```OpenBabel```, stored in ```output/receptor_files/```.

```Autodock Vina``` was executed using blind docking protocol with the following parameters:

- **Search space**: Defined by receptor center of mass and bounding box dimensions
- **Exhaustiveness**: 32 (standard computational power)
- **Number of modes**: 8 (conformational poses per docking run)
- **CPU threads**: -1 (use all available cores)
- **Random seed**: 1234 (for reproducibility)
- **Pose output**: All 8 poses saved for analysis

Each PHB oligomer (1-12 units) was docked against all 8 proteins, resulting in 96 independent docking runs and 768 total poses analyzed. Binding affinity (kcal/mol), RMSD lower bound (RMSD_LB), and RMSD upper bound (RMSD_UB) were extracted from Vina output files for all poses. Results were organized into a structured dataset for downstream statistical analysis, stored in ```output/docking_results_summary.csv```.

## Activity 3: Statistical Interpretation of Docking Results

**Notebook**: ```02_docking_results_analysis.ipynb```

### Objective
Perform comprehensive statistical analysis of the 704 docking poses to identify optimal PHB chain length, evaluate binding consistency, and assess pose quality using ```scipy```. The results were compared using PGA oligomers (1-12 units) as validation.

### Methodology

**1. Descriptive Statistics**

Mean, standard deviation, quartiles, skewness, and kurtosis were calculated for binding affinity distributions across all poses and grouped by protein type and PHB chain length.

**2. Affinity Classification**

Poses were classified into four binding strength categories based on established thresholds:

- **Strong**: < -8.0 kcal/mol
- **Moderate**: -8.0 to -6.0 kcal/mol
- **Weak**: -6.0 to -4.0 kcal/mol
- **Very Weak**: > -4.0 kcal/mol

**3. Chain Length Optimization**

One-way ANOVA was performed to test for significant differences in mean affinity across all 12 PHB chain lengths. Post-hoc pairwise comparisons were conducted using independent samples t-tests to identify the optimal chain length (highest binding affinity).

**4. Consistency Analysis**

Coefficient of variation (CV) was calculated for each protein: CV = (std / |mean|) × 100. Lower CV values indicate greater consistency in binding across all chain lengths and poses.

**5. Pose Quality Assessment**

RMSD spread was calculated as: RMSD_spread = RMSD_UB − RMSD_LB. This metric quantifies the structural variability of binding poses. Poses with spread < 5 Å are considered well-defined.

**6. Effect Size Calculation**

Cohen's d was computed to quantify the magnitude of differences between chain length groups, independent of sample size.

**7. Visualization**

Four-panel figure displaying: (a) distribution histogram of binding affinities, (b) box plots of affinity by protein, (c) box plots of affinity by PHB chain length, and (d) scatter plot of affinity vs. RMSD spread with linear regression.

## Activity 4: Mechanical and Structural Properties of PHA Polymer

**Notebook**: ```03_mechanical_properties.ipynb```

### Objective
Calculate and compare structural and mechanical properties of PHB and other biocompatible polymers (PGA, PLLA, P4HB, PHV) across oligomer chain lengths (1-12 units) to identify optimal chain length for medical-grade suture applications.

### Methodology

Five different polyester oligomers were synthesized in silico using ```RDKit```:

- **PHB** (Polyhydroxybutyrate): Monomer `C[C@H](CC(=O)O)O`; monomer length = 3.5 Å (PHA)
- **PHV** (Poly-3-hydroxyvalerate): Monomer `CC(CC(=O)O)O`; monomer length = 3.7 Å (PHA)
- **PGA** (Polyglycolic acid): Monomer `C(C(=O)O)O`; monomer length = 2.8 Å (standard suture material)
- **PLLA** (Poly-L-lactic acid): Monomer `C[C@H](C(=O)O)O`; monomer length = 3.0 Å (standard suture material)
- **P4HB** (Poly-4-hydroxybutyrate): Monomer `C(CC(=O)O)CO`; monomer length = 4.2 Å (PHA; flexible alternative)

Oligomers (n = 1 to 12 units) were constructed by iteratively applying the esterification reaction SMARTS pattern to the monomer. Each oligomer underwent 3D coordinate generation using ETKDG algorithm followed by MMFF94 energy minimization.

For each oligomer, 15 metrics were computed to characterize structural and mechanical properties:

- **Molecular Weight (MW)**: Total mass calculated using ```RDKit``` descriptors
- **Radius of Gyration (Rg)**: Average spatial extent from center of mass using `Descriptors3D.RadiusOfGyration()`
- **End-to-End Distance (R_ee)**: Distance between first and last heavy atoms extracted from 3D coordinates
- **Asphericity**: Shape descriptor (0 = spherical, 1 = rod-like) using `CalcAsphericity()`
- **LogP**: Lipophilicity calculated via Wildman-Crippen method (`Crippen.MolLogP()`)
- **Topological Polar Surface Area (TPSA)**: Polar surface area indicator
- **Ester Count & Density**: Number of ester bonds and bonds per 1000 Da (degradability proxy)
- **Rotatable Bonds**: Freely rotating bonds indicating flexibility
- **Chain Flexibility**: Ratio of rotatable bonds to heavy atoms
- **Solvent-Accessible Surface Area (SASA)**: Surface area accessible to solvent molecules
- **Persistence Length**: Characteristic length of chain rigidity = (Rg² / contour_length)
- **Methyl Branches**: Number of methyl side groups via SMARTS pattern matching
- **Shape Ratio**: Elongation measure = R_ee / Rg

Radius of gyration scaling was analyzed using power law fitting: **Rg = b × n^ν**

where ν is the scaling exponent indicating chain behavior:
- ν = 0.5 (ideal chain, theta solvent)
- ν = 0.6 (self-avoiding walk, good solvent)
- ν ≥ 0.7 (extended/semi-rigid chain)

Linear regression on log-transformed data (ln Rg vs. ln n) extracted slope (ν), intercept (b), and R² fit quality for each polymer.

Literature-based molecular weight thresholds were applied [^1]:

- **Low MW Brittle Limit**: < 9 kDa
- **Medical Grade Minimum**: 100 kDa (clinically acceptable)
- **High MW Suture Grade**: 240 kDa (optimal strength)

For each polymer, required monomer count was calculated: n_medical = MW_target / MW_monomer, with predicted Rg extrapolated via scaling law.

Properties were plotted across all five polymers:

1. Radius of Gyration (Rg) vs. chain length
2. Asphericity (linearity) vs. chain length
3. LogP (hydrophobicity) vs. chain length
4. Ester Density (degradability) vs. chain length

Scaling exponents were used to extrapolate Rg values to medical-grade molecular weights on double-logarithmic plots, with shaded zones demarcating brittle (<9 kDa) and medical-grade (100-240 kDa) regions.

-----
[^1]: Effect of Chain Stereoconfiguration on Poly(3-hydroxybutyrate) Crystallization Kinetics: https://pubs.acs.org/doi/pdf/10.1021/acs.biomac.2c00682

## END