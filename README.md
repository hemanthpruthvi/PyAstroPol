# PyAstroPol
Instrumental Polarization Analysis of Astronomical Optics

## Overview
The package has one simple goal : compute 4x4 Mueller matrix for the given optical system, and it is developed keeping astronomical optics in view.
It uses geometric optics approach i.e., all the analysis uses strictly ray treatment. Users should keep in mind that this is NOT for optical design i.e., it is presumed that the user already knows the optical system that is to be analyzed.

The package imports following external libraries, all of which are ubiquitous. They accompany any decent scientific `Python` distribution hence this package should function with virtually no dependency issues **however** it should be noted that it is developed with `Python3.6`.
```python
numpy
matplotlib
```

## Getting Started

### Installation

Follow these steps to start using the package.

1. Go to the user directory where the package is to be installed.  
`cd <User directory>`   
Download the package from the Github.   
`git clone https://github.com/hemanthpruthvi/PyAstroPol.git`  
Rename the top directory from `PyAstroPol.git` to `PyAstroPol`

2. Go to `PyAstroPol` root directory
`cd <User directory/PyAstroPol>`
Install the required dependencies by running  
`pip install -r requirements.txt`

3. Add `<User directory>/PyAstroPol` to the `PYTHONPATH` environment variable.  
In Windows systems, this option can be found at `Control Panel > All Control Panel Items > System > Advanced system settings > Environment Variables`   
In Linux systems, this can be done with the command line  
`export PYTHONPATH=<User directory>/PyAstroPol`

4. Import this package to your `Python` script using   
`import PyAstroPol as pap`

### Examples

[PyAstroPol/Examples/](https://github.com/hemanthpruthvi/PyAstroPol/tree/master/Examples) contains several examples files to demonstrate the applications of the package. They are provided in the form of `IPython` notebooks, and running them is a good way to quick-start using the package. **They also function as the test cases**. 

### Analysis on own

As previously mentioned, this is not a design software. Hence, one needs to know the optical system they wish to analyze. As per the framework of the `PyAstroPol` optical system, there are thee types of objects :
1. Source  
2. Components  
3. Detector   

Following steps illustrate how to devise a simple optical system.  
1. import the required modules.
```python
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import PyAstroPol as pap
```  
2. Create a source, and optionally create a source for display. For analysis one can define a source with a lot of rays (say 10000), and for display one can define a source with fewer rays (say 10).  
```python
S_analysis = pap.Source(10000, Clear=20)        # Source for analysis, with 10k rays and 20 mm size
S_display = pap.Source(10, Clear=20)            # Source for disply, with 10 rays and 20 mm size
```  
3. Create a component such as surface, lens etc., and position it. 
```python
L = pap.UncoatedLens(50, Thick=10, R1=200, R2=-200)     # Simple bi-convex lens of 50 mm size
L.translateOrigin(z=100.0)                              # Move the lens from default position (origin)
```  
4. Create a detector and position it.
```python
D = pap.Detector(50)                    # Detector of size 50 mm
D.translateOrigin(z=200.0)              # Move the detector from default position (origin)
```  
5. Put them together to create the optical system.
```python
O_system = pap.System(S_analysis, [L], D, dRays=S_display)
O_system.propagateRays()                                # Propagate rays in the optical system
```  
6. Display the optical system using matplotlib 3d axis.
```python
Fig = plt.figure()
Ax = Fig.add_subplot(111, projection='3d')
O_system.draw(Ax)
plt.show()
```  
7. Compute the Mueller matrix and print it.
```python
MM, T = O_system.getSystemMuellerMatrix()       # Compute Mueller matrix for the system
print(MM)
```

## Directories

[PyAstroPol/PyAstroPol/](https://github.com/hemanthpruthvi/PyAstroPol/tree/master/PyAstroPol)  
It is the main directory containing all the source files. **For more information on the `Classes` check the documentation hosted at [ReadTheDocs](https://pyastropol.readthedocs.io/en/latest/)**.

[PyAstroPol/Materials/](https://github.com/hemanthpruthvi/PyAstroPol/tree/master/Materials)  
It has the refractive index data for different materials in a formatted manner. These files are loaded by the code to look-up the refractive index information of the given material. Users can easily create such files using following steps.
1. Download wavelength vs refractive index file as `.csv` from popular refractive index database [RefractiveIndexInfo](https://refractiveindex.info/).
2. Rename the file to an appropriate material name e.g., for Aluminium the file name is `Al.csv`.
3. Copy the `.csv` file into `PyAstroPol/Materials/` directory.
4. Format the material file using provided function i.e., `formatMaterialFile(MaterialName)` (without file extensions).
5. The material is ready to be used by the code e.g., `M1 = Surface(50, n2='Al', Mirror=True)`.
6. An [example](https://github.com/hemanthpruthvi/PyAstroPol/blob/master/Examples/09_FormatMaterialFile.ipynb) is also provided.

[PyAstroPol/Docs/](https://github.com/hemanthpruthvi/PyAstroPol/tree/master/Docs)  
It contains documentation related codes and files.
`Theory_and_Implementation_Notes.ipynb` details the formulation behind the codes. **Users interested in development are encouraged to refer this document**.

## Conventions used in this package  
The most important aspect to remember while using the package is **the convention**, which is described below. 
### For astronomy : 
Positive X-axis : West  
Positive Y-axis : Zenith  
Positive Z-axis : North  
Positive Latitude : North  
Positive Hour Angle : West  
Positive Declination : North  
### For optics : 
Complex refractive index is **n-*i*k** where **n** and **k** are positive real numbers.    
Jones vector corresponding to positive Stokes-V is <img src="https://render.githubusercontent.com/render/math?math=\frac{1}{\sqrt 2} \begin{bmatrix} 1 \\ -i \end{bmatrix}">.

## Contributing
Any mode of contribution is highly encouraged.
1. Bug reporting : Open an issue in github with the following details.
    - Description of the bug
    - Python, numpy and matplotlib versions
    - Operating system details
    - Snippet of the code causing the issue
2. Feature request : Open an issue in github with the following details.
    - Description of the feature
    - Description of the application
    - If possible, an example
3. Example request : Open an issue in github with following details.
    - Description of the application
    - Expected output from the example
4. Bug fixes : Open a pull request in github with following details.
    - Description of the bug corresponding to the fix
5. Feature addition : Open a pull request in github with following details.
    - Description of the feature
    - Description of the application
    - At least one Example using the particular feature
6. Other : Open an issue on github with a description.

Kindly use appropriate Tags as well.

## TODO
1. Add polarizing elements such as birefringent elements, waveplates and polarizers.
2. Add feature to create and save coatings as files.
3. Add rectandular and elliptical apertures.