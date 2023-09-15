## molecular_fractions
A python script based on molecular formula and molecular number to calculate mass fraction, mole fraction...

### Usage

- calculating mass fraction:

  ```bash
  python cal_fraction.py -mass {A string containing the molecular formula and the number of molecules}
  ```

  - for example:

  ```bash
  python cal_fraction.py -mass "H2O 360 NaCl 4"
  ```

  - output on screen:

  ```bash
  Command line:
                  $python cal_fraction.py -mass H2O 360 NaCl 4
  
  # -------------------- Claculating mass fraction -------------------- #
  
  {'H2O': 360, 'NaCl': 4}
  
  Total mass of molecules = 6719.27188
  
  Mass fraction of H2O = 0.965209
  
  Mass fraction of NaCl = 0.034791
  
  # -------------------- Mass fraction end!!!!!!!! -------------------- #
  ```

  

- calculating mole fraction:

  ```bash
  python cal_fraction.py -mole {A string containing the molecular formula/symbol and the number of molecules}
  ```

  - for example:

  ```bash
  python cal_fraction.py -mole "H2O 368 CH4 64"
  ```

  - output:

  ```bash
  Command line:
                  $python cal_fraction.py -mole H2O 368 CH4 64
  
  # -------------------- Calculating mole fraction -------------------- #
  
  Mole fraction of H2O = 0.851852
  Mole fraction of CH4 = 0.148148
  
  # -------------------- Mole fraction end!!!!!!!! -------------------- #
  ```

  

