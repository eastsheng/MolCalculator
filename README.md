## molecular_fractions
A python script based on molecular formula and molecular number to calculate mass fraction, mole fraction...

### Usage

- calculating mass fraction:

  ```bash
  python cal_fraction.py -mass {A string containing the molecular formula and the number of molecules}
  ```

  - for example:

  ```bash
  python cal_fraction.py -mass "H2O 368 CH4 64"
  ```

  - output on screen:

  ```bash
  Command line:
                  $python cal_fraction.py -mass H2O 368 CH4 64
  
  # -------------------- Claculating mass fraction -------------------- #
  
  {'H2O': 368, 'CH4': 64}
  
  Total mass of molecules = 7656.340480000001
  
  Mass fraction of H2O = 0.8659
  
  Mass fraction of CH4 = 0.1341
  
  # -------------------- Mass fraction end!!!!!!!! -------------------- #
  ```

  

- calculating mole fraction:

  ```bash
  python cal_fraction.py -mole {A string containing the molecular formula/symbol and the number of molecules}
  ```

  - for example:

  ```bash
  python cal_fraction.py -mole "H2O 360 NaCl 4"
  ```

  - output:

  ```bash
  Command line:
                  $python cal_fraction.py -mole H2O 360 NaCl 4
  
  # -------------------- Calculating mole fraction -------------------- #
  
  H2O 摩尔分数: 0.989011
  NaCl 摩尔分数: 0.010989
  
  # -------------------- Mole fraction end!!!!!!!! -------------------- #
  ```

  
