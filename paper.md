---
title: '``idvf``: Iterative Inversion of Deformation Vector Field with Adaptive Bi-residual Feedback Control'
tags:
  - deformation vector field
  - optical flow
  - deformable registration
  - fixed-point iteration
  - inverse consistency
  - adaptive feedback control
  - medical image analysis
  - computer vision
  - MATLAB
authors:
  - name: Alexandros-Stavros Iliopoulos
    orcid: 0000-0002-1959-9792
    affiliation: 1
  - name: Abhishek Dubey
    orcid: 0000-0001-8052-7416
    affiliation: 1
  - name: Xiaobai Sun
    affiliation: 1
affiliations:
  - name: Department of Computer Science, Duke University, Durham, NC 27708, USA
    index: 1
date: 1 November 2018
bibliography: references.bib
---


# Summary

We provide a package for fast and accurate inversion of a deformation
vector field (DVF).  A DVF, also known as flow or dense motion field,
describes a non-linear mapping between a study image and a reference image.
DVFs are fundamental to several medical image analysis and computer vision
applications.  The inverse DVF is often needed together with the forward
DVF to enable back-and-forth mappings, and composite mappings among
multiple images.  Relevant applications include 4D image reconstruction,
anatomical atlas generation, dose accumulation estimation in adaptive
radiotherapy, simultaneous deformable registration, and symmetric
registration completion.  To our knowledge, no other package is available
for DVF inversion with guaranteed convergence.

We use an iterative inversion procedure and employ adaptive bi-residual
feedback control to achieve global convergence and local acceleration
[@medphys2018;phddubey2018].  The inverse DVF estimate and the forward DVF
must meet the inverse consistency (IC) condition.  We define two IC
residuals, which measure the inconsistency between two DVFs in the study
and reference domains.  The residuals, which are computationally available,
are related to the unknown inversion error.  We use them as feedback in a
two-phase iteration.  In the first phase, we modulate the study-domain IC
residual with adaptive feedback control for guaranteed global convergence,
under certain mild condition [@medphys2018].  Once the error is made
sufficiently small we switch to phase two, where we use the
reference-domain IC residual and achieve locally quadratic convergence
rate.  Phase transition and integration is enabled by a multi-resolution
scheme [@phddubey2018].

A pre-release version of ``idvf`` has been used in scientific publications
to demonstrate the inversion algorithm framework and analysis
[@medphys2018] and potential applications [@phddubey2018].  The ``idvf``
package is implemented in MATLAB.  Its source code is archived with Zenodo
[@idvfzenodo].



# References

