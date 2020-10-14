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

Polarization analysis is one of the key aspects of the optical system analysis in astronomy. The majority of the optical surfaces in the astronomical instruments induce a change in the state of polarization of the incoming light, which could lead to inaccurate measurements of the state of polarization or _polarimetry_. Hence, polarimetric instruments inevitably require the information of polarization properties of the optical system or the _instrumental polarization_. It can be determined through experiments, analysis or a combination of both. As the telescopes become larger and instruments become more complex, the instrumental polarization analysis has become all the more crucial (e.g.,  @TMT_2015 ; @DKIST_2016).  

# Summary

The Python package `PyAstroPol` provides a means to analyze the polarization properties of a given optical system with relative ease and minimal dependencies. The one simple goal is to calculate the Mueller matrix (e.g., @Gil_2016_MullMatAppr) of the optical system.  
It can be easily distributed along with the polarization calibration routines. Albeit there are several limitations, the available features should cover a large set of the polarization analysis related applications. The limitations and salient features are listed below.

Important Limitations :   
- All the analysis uses strictly the ray treatment. Hence, all the limitations of the rays optics shall be applicable.   
- Only circular optics (apertures) can be devised at the moment.   
- Incidentally, birefringent components are not included yet. The justification for this choice is that the behaviour of the birefringent components is fairly straightforward as they strongly polarize the light.   
- Limited visualization features.
- It is neither a design software nor interactive. That is, the user must know the optical system that is to be analyzed and the system must be updated every time a component is changed.

Salient features :   
- Calculation of the Mueller matrix of a given optical system is a simple and only end goal. All the rest are by-products of such analysis. The system Mueller matrix is computed by coherently adding the electric field vectors after propagation. It implies that the system is imaging type, which is commonplace in astronomy.   
- Astronomical sources can be directly placed using relevant coordinates, namely, declination, hour angle and latitude of the telescope site.   
- Off-axis components are facilitated, as they have a significant effect on polarization.   
- Effect of multi-layered coatings, such as oxide layers and protective coatings, on the state of polarization, is included. These are also significant in the polarization analysis of the astronomical optics (e.g., @VanHarten_2009).   
- All the data, such as points of incidence, polarization directions, complex electric field values and more, are readily available to the user for any further analysis.   
- Material refractive index information can be downloaded from the popular online source [https://refractiveindex.info/](https://refractiveindex.info/) as `.csv`. It can be formatted and used with this package.   
- Spot diagram is possible at any instance, as a by-product of ray tracing.   

Usually, such analyses are carried out by devising telescope specific polarization models using scientific programming languages (open-source/commercial) and/or by utilizing commercially available analysis software (e.g., @TMT_2015, @DKIST_2016). Here, I attempt to draw the best of the both: open-source tools developed using scientific programming language, with a wide range of applications in modelling the astronomical optics. To my knowledge, there is no open-source software available with such polarization analysis features.

Results of the code have been verified against previous works (@Pruthvi2018), and Zemax OpticStudio&reg; – a popular commercial software. A variety of examples have been provided, and they should facilitate the quick-start. One of the examples also provides the aforementioned comparison with the Zemax OpticStudio&reg;.

# Acknowledgements

This work has been carried out as a part of the ongoing project _Jets in the solar atmosphere_, funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) - Projektnummer 407727365. I thank DVS Phanindra (IIA, Bengaluru) and V. Sreekanth Reddy (CHESS, Hyderabad) for the discussion, and Mathias Waidele (KIS, Freiburg) for his feedback. I also thank the PI of the DFG project Markus Roth (KIS, Freiburg) for the opportunity to carry out this work. 

# References