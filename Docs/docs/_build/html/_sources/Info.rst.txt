About the Package
=================

Introduction
------------
The package ``PyAstroPol`` provides a simple way to perform instrumental polarization analysis of the astronomical optics. The package has one simple goal: compute 4x4 Mueller matrix for the given optical system. It is developed keeping astronomical optics in view. The geomatric optics approach is used in the package hence all the analysis uses strictly ray treatment. Users should keep in mind that the software is NOT for the optical design purposes i.e., it is presumed that the user already knows the optical system that is to be analyzed.

Keywords
--------
|   Python
|   Astronomy
|   Astronomical optics
|   Telescopes
|   Instrumentation
|   Polarization
|   Instrumental polarization
|   Polarimetry
|   Polarized ray tracing
|   Mueller matrix

Installation
------------
The package has minimal dependencies that accompany any decent scientific Python distribution. Hence, I believe that the package should function with virtually no dependency issues. However, it should be noted that it is developed with ``Python3.6``. Install these packages using the provided file ``requirements.txt``.
Follow these steps to start using the package.

|   1. Go to the selected user directory, with the coomand line ``cd <User directory>``
|      Download the package using ``git clone https://github.com/hemanthpruthvi/PyAstroPol.git``
|      Rename the directory from ``PyAstroPol.git`` to ``PyAstroPol``
|       
|      Alternatively, download the package as ``zip`` from ``https://github.com/hemanthpruthvi/PyAstroPol``
|      Extract it to the selected user directory directory as ``<User directory>/PyAstroPol``
|
|   2. Add ``<User directory>/PyAstroPol`` to the ``PYTHONPATH`` environronment variable.
|      In Windows systems, this option can be found at ``Control Panel > All Control Panel Items > System > Advanced system settings > Environment Variables``
|      In Linux systems, this can be done with the command line ``export PYTHONPATH=<User directory>/PyAstroPol``
|
|   3. Import this package to your ``Python`` script using ``import PyAstroPol``

Conventions
-----------
The most important aspect to remember while using the package is the convention, which is described below. 

Astronomy
    | Positive X-axis : West
    | Positive Y-axis : Zenith
    | Positive Z-axis : North 
    | Positive Latitude : North 
    | Positive Hour Angle : West 
    | Positive Declination : North 

Optics
    | Complex refractive index is **n**-*i*.**k** where **n** and **k** are positive real numbers.    
    | Jones vector corresponding to positive Stokes-V is :math:`\frac{1}{\sqrt{2}} \begin{bmatrix} 1 \\ -i \end{bmatrix}`
