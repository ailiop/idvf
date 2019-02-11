# idvf: _Iterative inversion of deformation vector field with adaptive bi-residual feedback control_


[![Zenodo DOI](https://zenodo.org/badge/153499710.svg)](https://zenodo.org/badge/latestdoi/153499710)
[![GitHub release](https://img.shields.io/github/release/ailiop/idvf.svg)](https://github.com/ailiop/idvf/releases/)
[![GitHub license](https://img.shields.io/github/license/ailiop/idvf.svg)](https://github.com/ailiop/idvf/blob/master/LICENSE)
[![Github all releases](https://img.shields.io/github/downloads/ailiop/idvf/total.svg)](https://github.com/ailiop/idvf/releases/)
[![GitHub issues](https://img.shields.io/github/issues/ailiop/idvf.svg)](https://github.com/ailiop/idvf/issues/)




<a name="contents"></a>

## Contents


- [Algorithm overview](#algorithm)
	- [Two inverse consistency residuals](#ic-residuals)
	- [DVF inversion iteration with bi-residual feedback control](#inversion-iteration)
	- [References](#references)
- [Getting started](#getting-started)
	- [Installation](#installation)
	- [Testing](#testing)
- [Software description](#software)
	- [DVF inversion](#inversion-function)
	- [Preprocessing & DVF characterization](#preprocessing)
	- [Post-evaluation](#post-evaluation)
	- [GPU utilization](#gpu)
	- [Eigenvalue calculations](#eigenvalues)
- [System environment](#system-reqs)
- [License and community guidelines](#license-contrib-reports)
- [Contributors](#contributors)




<a name="algorithm"></a>

## Algorithm overview


This repository contains a set of MATLAB functions for fast and accurate
inversion of 2D/3D deformation vector fields (DVFs), which are also known
as flow or dense motion fields.  The inversion algorithm framework and a
unified analysis are described in References
[[1](#medphys2018),[2](#phddubey2018)].  To our knowledge, this is the
first DVF inversion framework that is supported by theoretical guarantees
under mild and verifiable conditions.  We give a brief overview of the DVF
inversion algorithm which underlies the repository functions.


We assume that we are provided with a deformation vector field (DVF)
denoted by _u_.  The DVF describes a (non-linear) mapping from a reference
image _I<sub>r</sub>_ onto a study image _I<sub>s</sub>_ via point-wise
displacement, i.e.,

<p align="center"><i>
	I<sub>r</sub>(y) = I<sub>s</sub>(y + u(y)),
</i></p>

at all locations _y_ in the reference domain.  We aim at computing the
inverse DVF _v_ such that

<p align="center"><i>
	I<sub>s</sub>(x) = I<sub>r</sub>(x + v(x)),
</i></p>

at all _x_ in the study domain.



<a name="ic-residuals"></a>

### Two inverse consistency residuals


If two DVFs _u_ and _v_ are inverse of each other, they must meet the
*inverse consistency* (IC) condition.  We define the study IC residual,

<p align="center"><i>
	s(x) = v(x) + u(x + v(x)),
</i></p>

and the reference IC residual,

<p align="center"><i>
	r(y) = u(y) + v(y + u(y)).
</i></p>

The DVFs _u_ and _v_ are inverse-consistent if both residuals are zero.
The residuals are related to the inversion error, _e(x) = v(x) &minus;
v&ast;(x)_, where _v&ast;_ is the true inverse DVF.  The IC residuals are
computationally available, while the inversion error is unknown.  We use
the IC residuals as feedback in iterative DVF inversion, and exercise
adaptive control over the feedback for global convergence and local
acceleration.



<a name="inversion-iteration"></a>

### DVF inversion iteration with bi-residual feedback control


Due to the non-linear nature of the DVF mapping, one resorts to iterative
procedures for DVF inversion.  Our procedure has two stages: preprocessing
and adaptive inversion iteration.


In the preprocessing stage, we compute the DVF Jacobians and their
eigenvalues over the spatial domain.  We use this spectral information to
determine data-adaptive parameters for feedback control.


The inversion iteration proceeds in two phases.  During phase one, at step
_k+1_, we update the current inverse estimate _v<sub>k</sub>_ using the
study IC residual _s<sub>k</sub>_ as feedback:

<p align="center"><i>
	v<sub>k+1</sub>(x) = 
	v<sub>k</sub>(x) &minus; (1 &minus; μ<sub>k</sub>(x)) s<sub>k</sub>(x),
</i></p> 

where _μ_ is the feedback control parameter.  With the adaptive control
schemes provided in this repository, the iteration is guaranteed to
converge globally, under certain mild conditions [[1](#medphys2018)].


Once the errors are made sufficiently small, the iteration is switched to
phase two:

<p align="center"><i>
	v<sub>k+1</sub>(x) = 
	v<sub>k</sub>(x) &minus; r<sub>k</sub>(x + v<sub>k</sub>(x)), 
</i></p>

where the reference IC residual _r<sub>k</sub>_ at displaced locations is
used as feedback.  Inversion errors can be estimated using the IC residuals
and the DVF spectral information.  The local convergence rate is quadratic
during phase two.  We may refer to this iteration as implicit Newton
iteration.  However, phase-two iteration steps do not entail explicit
formation and inversion of the Newton matrix; the explicit Newton step
suffers from several numerical problems.


Phase transition and integration in the inversion iteration is facilitated
by a multi-resolution scheme, elaborated in Reference [[2](#phddubey2018)].



<a name="references"></a>

### References


<a name="medphys2018"></a>

[1] A. Dubey&ast;, A.S. Iliopoulos&ast;, X. Sun, F.F. Yin, and L. Ren,
["Iterative inversion of deformation vector fields with feedback
control"][medphys2018-doi], *Medical Physics*, 45(7):3147-3160, 2018.
<small>[(arXiv)][medphys2018-arxiv]</small>

[medphys2018-doi]:   https://doi.org/10.1002/mp.12962

[medphys2018-arxiv]: https://arxiv.org/abs/1610.08589


<a name="phddubey2018"></a>

[2] A. Dubey, *Symmetric completion of deformable registration via
bi-residual inversion*, PhD thesis, Duke University, Durham, NC, USA, 2018.




<a name="getting-started"></a>

## Getting started



<a name="installation"></a>

### Installation


To use idvf, simply add its top-level directory to the [MATLAB
path][matlab-path].  All functions are organized in
[packages][matlab-packages].

[matlab-path]:     https://www.mathworks.com/help/matlab/matlab_env/what-is-the-matlab-search-path.html

[matlab-packages]: https://www.mathworks.com/help/matlab/matlab_oop/scoping-classes-with-packages.html


<a name="testing"></a>

### Testing


To test idvf, run the included demo scripts (`demo_inversion_2d`,
`demo_inversion_3d_z0`, `demo_inversion_3d_zsin`).  Each script
demonstrates inversion of a sample forward DVF using 8 different feedback
control settings, and produces the following visualizations:

-   Deformation of a synthetic grid-like image by the forward DVF.
-   Spectral measure maps of the forward DVF
    ([preprocessing](#preprocessing)).
-   IC residual magnitudes ([post-evaluation](#post-evaluation)) with each
    control setting: percentiles at each step of the iteration, and spatial
    maps at the end of the iteration.
-   Image-space difference maps in reference image recovery with each
    control setting.


The adaptive control schemes are guaranteed to converge; that is, the
corresponding IC residuals tend towards zero, at a rate that depends on the
control setting.  The spectral measure maps provide useful information
regarding the invertibility and controllability of the forward DVF over its
spatial domain.




<a name="software"></a>

## Software description



<a name="inversion-function"></a>

### DVF inversion


The main function is `dvf.inversion`.  It takes a forward DVF _u_ and
returns the inverse DVF _v_.  A 3D DVF is represented by a 4D array of size
_N<sub>x</sub>_ x _N<sub>y</sub>_ x _N<sub>z</sub>_ x _3_, in `single` or
`double` precision.  Specifically, `U(i,j,k,1:3)` is the 3D forward
displacement vector _u(x)_ at voxel _x = (i,j,k)_.  Two-dimensional DVFs
are represented analogously.


We provide a number of options in function `dvf.inversion`, which are input
as name-value pairs.  These options pertain to the choice of feedback
control schemes and iteration parameters.  The function documentation
details the supported options.


The memory requirement of `dvf.inversion` is approximately 5 times the
amount of memory for the input DVF data array.



<a name="preprocessing"></a>

### Preprocessing & DVF characterization


Preprocessing functions are invoked automatically by the main function when
adaptive feedback control is chosen.  Functions `dvf.jacobian` and
`dvf.eigJacobian` compute the DVF Jacobians and their eigenvalues,
respectively, at all pixels/voxels.  Function `dvf.feedbackControlVal`
takes the eigenvalues and returns adaptive feedback control parameter
values per the chosen scheme.


Function `dvf.ntdcMeasures` calculates three characteristic spectral
measures of a DVF.  Specifically, it returns the algebraic control index
map, the non-translational component spectral radius map, and the
determinant map.



<a name="post-evaluation"></a>

### Post-evaluation


The two IC residuals enable post-evaluation of inverse DVF estimates in the
absence of the true inverse DVF.  Function `dvf.inversion` returns IC
residual magnitude percentiles at each iteration step as optional output.
Function `dvf.icResidual` returns an IC residual field given a DVF pair
(reference or study IC residual, depending on the order of the input DVFs).


Certain regions of the study domain are invalid for IC residual
calculations.  These regions cover points that are either mapped outside
the domain by the forward DVF, or lie outside the deformed reference
domain.  Function `dvf.maskDomain` calculates the valid domain.



<a name="gpu"></a>

### GPU utilization


All `dvf.*` functions are GPU-enabled, supported by the MATLAB Parallel
Computing Toolbox.  When a compatible GPU is available, GPU computation is
invoked by casting the input data arrays into `gpuArray` objects; for
example:

    V = dvf.inversion(gpuArray(U));
	
The output data arrays are returned as `gpuArray` objects.



<a name="eigenvalues"></a>

### Eigenvalue calculations


Eigenvalue calculations are currently slow with large DVFs.  They are
implemented as a loop that calls the MATLAB `eig` function for each
Jacobian over the spatial domain.  This approach is inefficient but has the
benefit of numerical stability.  In our experience with respiratory DVFs
over thoracic and abdominal CT domains, we have not encountered any major
issues with efficient but numerically unstable alternatives such as
closed-form root calculation.  Nevertheless, we do not include such
alternatives, in the interest of stability and simplicity over efficiency.




<a name="system-reqs"></a>

## System environment


The repository code was developed and tested on MATLAB R2018b.  The
following functions in MATLAB toolboxes are used:

-   `imresize`, `imresize3`, `imdilate`, `imerode`, and `imclose` (Image
    Processing Toolbox, tested with version 10.3);
-   `prctile` (Statistics and Machine Learning Toolbox, tested with version
    11.4); and
-   `gpuArray` (Parallel Computing Toolbox, tested with version 6.13).


GPU computations were tested on an NVIDIA Quantum TXR113-1000R and CUDA
10.0 drivers.




<a name="license-contrib-reports"></a>

## License and community guidelines


The idvf code is licensed under the [GNU general public license
v3.0][license].  If you wish to contribute to idvf or report any
bugs/issues, please see our [contribution guidelines][contrib] and [code of
conduct][conduct].

[license]: https://github.com/ailiop/idvf/blob/master/LICENSE

[contrib]: https://github.com/ailiop/idvf/blob/master/CONTRIBUTING.md

[conduct]: https://github.com/ailiop/idvf/blob/master/CODE_OF_CONDUCT.md



<a name="contributors"></a>

## Contributors


-   *Design, development, testing:*  
    Alexandros-Stavros Iliopoulos, Abhishek Dubey, and Xiaobai Sun  
    Department of Computer Science, Duke University

-   *Consulting on medical CT image analysis and applications:*  
    Lei Ren and Fang-Fang Yin  
    Department of Radiation Oncology, Duke University Medical School

-   *Review, suggestions, bug reports:*  
	@Kevin-Mattheus-Moerman
