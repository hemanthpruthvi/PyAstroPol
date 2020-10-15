---
title: 'PyAstroPol : A Python package for the instrumental polarization analysis of the astronomical optics.'

tags:
  - Python
  - Astronomy
  - Polarization
  - Optics

authors:
  - name: Hemanth Pruthvi.
    orcid: 0000-0002-4892-6561
    affiliation: "1"

affiliations:
 - name: Leibniz-institut für Sonnenphysik, Freiburg, Germany.
   index: 1

date: 25 August 2020

bibliography: paper.bib
---

# Statement of Need

The instrumental polarization analysis is one of the key aspects of the optical system analysis of the astronomical telescopes and instruments. The majority of the optical surfaces in the astronomical instruments induce a change in the state of polarization of the incoming light, which could lead to inaccurate measurements of the state of polarization or _polarimetry_. Hence, polarimetric instruments inevitably require the information of polarization properties of the optical system or the _instrumental polarization_. It can be determined through experiments, analysis or a combination of both. As the telescopes become larger and instruments become more complex, the instrumental polarization analysis has become all the more crucial (e.g.,  @TMT_2015 ; @DKIST_2016).  

# Summary

The Python package `PyAstroPol` provides a means to analyze the polarization properties of a given optical system with relative ease and minimal dependencies. The one simple goal is to calculate the Mueller matrix (e.g., @Gil_2016) of the optical system.  
It can be easily distributed along with the polarization calibration routines. 

In the polarization analysis of the astronomical telescopes, various approaches have been adopted considering the complexity of the system. A significant part of the complexity is due to the fact that the instrumental polarization is often time-dependent. For the solar telescopes with Coelostats, Mueller matrices were analytically derived as a function of time (e.g., @KTT_1985, @VTT_2005). For Thirty Meter Telescope (TMT), a combination of analytical and numerical methods is used (@TMT_2015). However, for Daniel K. Inouye Solar Telescope (DKIST), Zemax – a commercial software is used, along with the team's in-house tools (@DKIST_2016).

`PyAstroPol` aims to provide better parts of the aforementioned approaches: open-source tools developed using scientific programming language, with a range of applications in modelling the astronomical optics. There is no open-source software which has polarization ray propagation features as per my knowledge. Hence, Zemax OpticStudio&reg; modelling is used for the comparison. The salient features and important limitations are listed below.

Salient features :   
- Calculation of the Mueller matrix of a given optical system is a simple and only end goal. Rest are the by-products of such analysis. The system Mueller matrix is computed by coherently adding the electric field vectors after propagation. It implies that the system is imaging type, which is commonplace in astronomy.   
- Astronomical sources can be directly placed in the model using relevant coordinates, namely, declination, hour angle and latitude of the telescope site.   
- Off-axis components are facilitated, as they have a significant effect on polarization.   
- Effect of multi-layered coatings, such as oxide layers and protective coatings, on the state of polarization, is included. These are also significant in the polarization analysis of the astronomical optics (e.g., @VanHarten_2009).   
- All the data, such as points of incidence, polarization directions, complex electric field values and more, are readily available to the user for any further analysis.   
- Material refractive index information can be downloaded from the popular online source [https://refractiveindex.info/](https://refractiveindex.info/) as `.csv`. It can be formatted and used with this software.   
- Spot diagram is possible at any instance, as a by-product of ray tracing.

Important Limitations :   
- All the analysis uses strictly the ray treatment. Hence, all the limitations of the rays optics shall be applicable.   
- Only circular optics (apertures) can be devised at the moment.   
- Incidentally, birefringent components are not included yet. The justification for this choice is that the behaviour of the birefringent components is fairly straightforward as they strongly polarize the light.   
- The visualization features are limited.
- It is neither a design software nor interactive. That is, the user must know the optical system that is to be analyzed, and the system must be updated every time a component is changed.

Results of the code have been verified against previous works (@Pruthvi_2018), and Zemax OpticStudio&reg;. A variety of examples have been provided, and they should facilitate the quick-start. One of the examples also provides the aforementioned comparison with the commercial software. 

# Acknowledgements

This work has been carried out as a part of the ongoing project _Jets in the solar atmosphere_, funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) – Projektnummer 407727365. I thank DVS Phanindra (IIA, Bengaluru) and V. Sreekanth Reddy (CHESS, Hyderabad) for the discussion, and Mathias Waidele (KIS, Freiburg) for his feedback. I also thank the PI of the DFG project Markus Roth (KIS, Freiburg) for the opportunity to carry out this work. 

# References