# PyAstroPol
Instrumental Polarization Analysis for Astronomical Optics

## Overview and Structure
The package has one simple goal : Compute 4x4 Mueller matrix for the given optical system, and it is developed keeping astronomical optics in view.
It uses geomatric optics approach i.e., all the analysis uses strictly ray treatment. The most important aspect to remember while delving into using the package is __the convention__, which is described in this file.

## Conventions used in the code

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

