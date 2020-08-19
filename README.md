# PyAstroPol
Instrumental Polarization Analysis for Astronomical Optics

## Overview
The package has one simple goal : Compute 4x4 Mueller matrix for the given optical system, and it is developed keeping astronomical optics in view.
It uses geomatric optics approach i.e., all the analysis uses strictly ray treatment. The most important aspect to remember while delving into using the package is __the convention__, which is described in this file.

## Directories
[./PyAstroPol/](https://github.com/hemanthpruthvi/PyAstroPol/tree/master/PyAstroPol)  
It is the main directory containing all the source files. For example, they can be imported as  
```python
from PyAstroPol.PyAstroPol import * 
R = Rays(10)
```
[./Materials/](https://github.com/hemanthpruthvi/PyAstroPol/tree/master/Materials)  
It has the refractive index data for different materials in a formatted manner. These files are loaded by the code to look-up the refractive index information of the given material. Users cas easily create such files using following easy steps.
1. Download wavelength vs refractive index file as `.csv` from popular refractive index database [RefractiveIndexInfo](https://refractiveindex.info/).
2. Rename the file to a material name which you can call e.g., for Aluminium the file name is `Al.csv`.
3. Copy the `.csv` file into `./Materials/` directory.
4. Format the material file using provided function i.e., `formatMaterialFile(MaterialName)` (without file exyensions).

[./Examples/](https://github.com/hemanthpruthvi/PyAstroPol/tree/master/Examples)  
It contains reveral examples files to demonstrate the applications of the package. They are provided in the form of `iPython` notebook and it is a good way to quick-start using the package.

[./Docs/](https://github.com/hemanthpruthvi/PyAstroPol/tree/master/Docs)  
It contains documentation related codes and files.  









## Conventions used in this package

### For astronomy : 
Positive X-axis : West  
Positive Y-axis : Zenith  
Positive Z-axis : North  

Positive Latitude : North  
Positive Hour Angle : West  
Positive Declination : North  

### For optics : 
Complex refractive index is __n-__*i*__k__ where __n__ and __k__ are positive real numbers.    
Jones vectors corresponding to positive Stokes-V are <img src="https://render.githubusercontent.com/render/math?math=\frac{1}{\sqrt 2} \begin{bmatrix} 1 \\ 0 \end{bmatrix}"> and <img src="https://render.githubusercontent.com/render/math?math=\frac{1}{\sqrt 2} \begin{bmatrix} 0 \\ -i \end{bmatrix}">.  

