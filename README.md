# PyFORKIT


This repository contains simple script for reading MigroMag files and plotting FORC (first-order reversal curves) diagram as well as related figures. This script was inspired by FORKIT program which was written  Gary Acton (http://paleomag.ucdavis.edu/software-forcit.html)


How to run
----------------

```
python pyforkit.py -i path_to_file.frc
```

Another option can be found by 

```
python pyforkit.py -h
```


Required libraries
-------------------

This programm uses NumPy, SciPy and Matplotlib. It was tested against NumPy v1.14.0,  SciPy  v1.0.0 and Matplotlib v2.1.2 for Windows 10 and MacOS X.