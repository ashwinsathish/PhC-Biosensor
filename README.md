# Numerical modelling of a highly sensitive surface plasmon sensor using silicon and platinum diselenide stacks

This repository contains the implementation of a one-dimensional photonic crystal (1D PhC) based Surface Plasmon Resonance (SPR) sensor designed for biomolecular detection. The sensor utilizes aluminum as the plasmonic metal, combined with silicon and platinum diselenide layers for enhanced sensing capabilities.

## Features

- High-sensitivity SPR sensing in the near-infrared region (1550 nm)
- Label-free detection of biomolecular interactions
- Optimized multi-layer structure using Al/Si/PtSe2
- Transfer matrix method implementation for optical simulations
- Angle interrogation analysis capabilities

## Sensor Architecture

The sensor structure consists of:
- Base: CaF2 prism
- Plasmonic layer: 30nm Aluminum
- Periodic structure: 3 stacks of
  - Silicon (40nm)
  - Platinum diselenide (4nm)
- Analyte layer: 2500nm

## Performance Metrics

- Sensitivity: 101.1°/RIU
- Figure of Merit (FOM): 1531.82 RIU⁻¹
- Operating wavelength: 1550nm
- Refractive index detection range: 1.33-1.34

## Optimization

The repository includes scripts for:
- Si layer thickness optimization (10-60nm)
- PtSe2 thickness optimization (2-10nm)
- Stack number optimization (0-6 stacks)
- Angle interrogation analysis

## References

1. **Panda, Abinash, Pukhrambam, Puspa, Dadure, Pankaj, & Pakray, Dr. Partha.** (2022). *Application of Machine Learning for Accurate Detection of Hemoglobin Concentrations Employing Defected 1D Photonic Crystal.* _Silicon, 14_. [https://doi.org/10.1007/s12633-022-01926-x](https://doi.org/10.1007/s12633-022-01926-x)  

2. **Shukla, S., Grover, N., & Arora, P.** (2023). *Resolution enhancement using a multi-layered aluminum-based plasmonic device for chikungunya virus detection.*  _Opt Quant Electron, 55_, 274. [https://doi.org/10.1007/s11082-022-04485-y](https://doi.org/10.1007/s11082-022-04485-y)  

3. **Uniyal, A., Pal, A., & Chauhan, B.** (2023). *Long-range SPR sensor employing platinum diselenide and cytop nanolayers giving improved performance.* _Physica B: Condensed Matter, 649_, 414487. [https://doi.org/10.1016/j.physb.2022.414487](https://doi.org/10.1016/j.physb.2022.414487)  

4. **Panda, A., & Pukhrambam, P. D.** (2022). *Modeling of High-Performance SPR Refractive Index Sensor Employing Novel 2D Materials for Detection of Malaria Pathogens.* _IEEE Transactions on Nanobioscience, 21_(2), 312-319. [https://doi.org/10.1109/TNB.2021.3115906](https://doi.org/10.1109/TNB.2021.3115906)  
