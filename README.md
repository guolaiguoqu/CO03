# Strange Particle Analysis (CDF Data)

This project analyses strange particle decays from proton–antiproton collision data (CDF experiment) 
to reconstruct particle masses and estimate lifetimes using kinematic information and selection criteria.

The analysis is implemented in MATLAB with a modular structure, separating data handling, reconstruction, and fitting. 
The reconstructed masses of $K^0_S$ and $\Lambda^0$ are consistent with known values within uncertainties, 
while lifetime estimates show larger systematic uncertainty.

The full methodology and results are presented in `report/report.pdf`.

## Structure

* `analysis/` – main data analysis and reconstruction
* `mathematical methods/` – track modelling, vertex reconstruction and visualisation
* `raw data wrapment/` – data handling utilities
* `report/` – report and results

## Notes

This work was developed as part of an undergraduate practical. 
Some data and supporting materials were provided; the analysis and code structure were implemented independently.

## Acknowledgements

Based on the practical *“Strange mesons from proton–antiproton collisions at the CDF experiment”*.
Further details are provided in "ACKNOWLEDGEMENTS.md".
