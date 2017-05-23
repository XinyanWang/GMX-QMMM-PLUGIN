# SIMPLE-QMMM
This project is made for using QMMM method in searching for reaction pathway on enzyme system.

To Do List

- reader
  - GroFile
  - TopFile
- engine (renew positions)
  - NVEEngine
    - constraint
  - LangevinEngine
  - NPTEngine
- calculator (calculate forces and velocities)
  - QMCalculator
    - DFTB+
    - Gaussian
  - MMCalculator
- reporter
  - PDBReporter
  - TRJReporter
- add force
  - ForceAdder
  - SimpleHarmonicAdder
  - GaussAdder
  - LambdaAdder
  - PrincipalComponentAdder
