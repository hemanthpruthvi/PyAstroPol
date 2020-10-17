About the Package
=================

The ``Python`` package ``PyAstroPol`` provides a simple way to perform instrumental polarization analysis of the astronomical optics. The package has one simple goal: compute 4x4 Mueller matrix for the given optical system. It is developed keeping astronomical optics in view. The geometric optics approach is used in the package hence all the analysis uses strictly ray treatment. Users should keep in mind that the software is NOT for the optical design purposes i.e., it is presumed that the user already knows the optical system that is to be analyzed. As it has minimal dependencies it can be easily used, or integrated/distributed with polarization calibration routines of the users.

``PyAstroPol`` aims to provide better aspects of the existing polarization analysis approaches: open-source tools developed using scientific programming language, with a range of applications in modelling the astronomical optics.  A variety of examples have been provided, and they should facilitate the quick-start. One of the examples also provides the comparison of the results with those of a popular commercial software Zemax OpticStudio®. 

Salient features
----------------

- Calculation of the Mueller matrix of a given optical system is a simple and only end goal. Rest are the by-products of such analysis. The system Mueller matrix is computed by coherently adding the electric field vectors after propagation. It implies that the system is imaging type, which is commonplace in astronomy.   
- Astronomical sources can be directly placed in the model using relevant coordinates, namely, declination, hour angle and latitude of the telescope site.   
- Off-axis components are facilitated, as they have a significant effect on polarization.   
- Effect of multi-layered coatings, such as oxide layers and protective coatings, on the state of polarization, is included. These are also significant in the polarization analysis of the astronomical optics.   
- All the data, such as points of incidence, polarization directions, complex electric field values and more, are readily available to the user for any further analysis.   
- Material refractive index information can be downloaded from the popular online source `RefractiveIndexInfo <https://refractiveindex.info/>`_ as ``.csv``. It can be formatted and used with this software.   
- Spot diagram is possible at any instance, as a by-product of ray tracing.

Important Limitations
---------------------

- All the analysis uses strictly the ray treatment. Hence, all the limitations of the rays optics shall be applicable.   
- Only circular optics (apertures) can be devised at the moment.   
- Incidentally, birefringent components are not included yet. The justification for this choice is that the behaviour of the birefringent components is fairly straightforward as they strongly polarize the light.   
- The visualization features are limited.
- It is neither a design software nor interactive. That is, the user must know the optical system that is to be analyzed, and the system must be updated every time a component is changed.


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

