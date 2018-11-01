# idvf

## _Iterative Inversion of Deformation Vector Field<br/>with Adaptive Bi-residual Feedback Control_


<a id="toc"></a>

### Table of contents

- [Algorithm overview](#algorithm)
  - [Two inverse consistency residuals](#ic-residuals)
  - [DVF inversion iteration with bi-residual feedback control](#inversion-iteration)
  - [References](#references)
- [Software description](#software)
  - [DVF inversion](#inversion-function)
  - [Preprocessing & DVF characterization](#preprocessing)
  - [Post-evaluation](#post-evaluation)
  - [GPU utilization](#gpu)
- [System environment](#system-reqs)
- [Contributors](#contributors)



<a id="algorithm"></a>

### Algorithm overview

This repository contains a set of MATLAB functions for fast and accurate
inversion of 2D/3D deformation vector fields (DVFs), which are also known
as flow or dense motion fields. The inversion algorithm framework and a
unified analysis are described in References [<a
href="#medphys2018">1</a>,<a href="#dukephd2018">2</a>]. We give a brief
overview of the DVF inversion algorithm which underlies the repository
functions.

We assume that we are provided with a deformation vector field (DVF)
denoted by _u_. The DVF describes a (non-linear) mapping from a reference
image _I<sub>r</sub>_ onto a study image _I<sub>s</sub>_ via point-wise
displacement, i.e.,

<p align="center"><i>
	I<sub>r</sub>(y) = I<sub>s</sub>(y + u(y)),
</i></p>

at all locations _y_ in the reference domain. We aim at computing the
inverse DVF _v_ such that

<p align="center"><i>
	I<sub>s</sub>(x) = I<sub>r</sub>(x + v(x)),
</i></p>

at all _x_ in the study domain.



<a id="ic-residuals"></a>

#### Two inverse consistency residuals

If two DVFs _u_ and _v_ are inverse of each other, they must meet the
*inverse consistency* (IC) condition. We define the study IC residual,

<p align="center"><i>
	s(x) = v(x) + u(x + v(x)),
</i></p>

and the reference IC residual,

<p align="center"><i>
	r(y) = u(y) + v(y + u(y)).
</i></p>

The DVFs _u_ and _v_ are inverse-consistent if both residuals are zero. The
residuals are related to the inversion error, _e(x) = v(x) - v&ast;(x)_,
where _v&ast;_ is the true inverse DVF. The IC residuals are
computationally available, while the inversion error is unknown. We use the
IC residuals as feedback in iterative DVF inversion, and exercise adaptive
control over the feedback for global convergence and local acceleration.



<a id="inversion-iteration"></a>

#### DVF inversion iteration with bi-residual feedback control

Due to the non-linear nature of the DVF mapping, one resorts to iterative
procedures for DVF inversion. Our procedure has two stages: preprocessing
and adaptive inversion iteration.

In the preprocessing stage, we compute the DVF Jacobians and their
eigenvalues over the spatial domain. We use this spectral information to
determine data-adaptive parameters for feedback control.

The inversion iteration proceeds in two phases. During phase one, at step
_k+1_, we update the current inverse estimate _v<sub>k</sub>_ using the
study IC residual _s<sub>k</sub>_ as feedback:

<p align="center"><i>
	v<sub>k+1</sub>(x) = 
	v<sub>k</sub>(x) &minus; (1 - μ<sub>k</sub>(x)) s<sub>k</sub>(x),
</i></p> 

where _μ_ is the feedback control parameter. With the adaptive control
schemes provided in this repository, the iteration is guaranteed to
converge globally, under certain mild conditions [<a
href="#medphys2018">1</a>].

Inversion errors can be estimated using the IC residuals and the DVF
spectral information. Once the errors are made sufficiently small, the
iteration is switched to phase two:

<p align="center"><i>
	v<sub>k+1</sub>(x) = 
	v<sub>k</sub>(x) &minus; r<sub>k</sub>(x + v<sub>k</sub>(x)), 
</i></p>

where the reference IC residual _r<sub>k</sub>_ at displaced locations is
used as feedback. The local convergence rate is quadratic during phase
two. We may refer to this iteration as implicit Newton iteration. However,
phase-two iteration steps do not entail explicit formation and inversion of
the Newton matrix; the explicit Newton step suffers from several numerical
problems.

Phase transition and integration in the inversion iteration is enabled by a
multi-resolution scheme, elaborated in Reference [<a
href="#dukephd2018">2</a>].



<a id="references"></a>

#### References

<a id="medphys2018"></a>[1] A. Dubey&ast;, A.S. Iliopoulos&ast;, X. Sun,
F.F. Yin, and L. Ren, <a
href="http://dx.doi.org/10.1002/mp.12962">Iterative inversion of
deformation vector fields with feedback control</a>, *Medical Physics*,
45(7):3147-3160 (2018).

<a id="dukephd2018"></a>[2] A. Dubey, *Symmetric Completion of Deformable
Registration via Bi-residual Inversion*, PhD dissertation, Duke University,
Durham, NC, USA, 2018.



<a id="software"></a>

### Software description



<a id="inversion-function"></a>

#### DVF inversion

The main function is `dvf.inversion`. It takes a forward DVF _u_ and
returns the inverse DVF _v_. A 3D DVF is represented by a 4D array of size
_N<sub>x</sub>_ x _N<sub>y</sub>_ x _N<sub>z</sub>_ x _3_, in `single` or
`double` precision. Specifically, `U(i,j,k,1:3)` is the 3D forward
displacement vector _u(x)_ at voxel _x = (i,j,k)_. Two-dimensional DVFs are
represented similarly.

We provide a number of options in function `dvf.inversion`, which are input
as name-value pairs. These options are for the choice of feedback control
schemes and iteration parameters. The function documentation details the
provided options.

The memory requirement of `dvf.inversion` is approximately 5 times the
amount of memory for the input DVF data array.

We provide two scripts, `demo_inversion_2d` and `demo_inversion_3d`, to
demonstrate how to use the repository functions. The scripts include
post-evaluation of the inversion results in terms of IC residual magnitude
percentiles throughout the iteration, per several control schemes.



<a id="preprocessing"></a>

#### Preprocessing & DVF characterization

Preprocessing functions are invoked automatically by the main function when
adaptive feedback control is chosen. Functions `dvf.jacobian` and
`dvf.eigJacobian` compute the DVF Jacobians and their eigenvalues,
respectively, at all pixels/voxels. The function `dvf.feedbackControlVal`
takes the eigenvalues and returns adaptive values for the feedback control
parameter _μ_ per the chosen scheme.

The function `dvf.ntdcMeasures` calculates three characteristic spectral
measures of a DVF. Specifically, it returns the algebraic control index
map, the non-translational component spectral radius map, and the
determinant map.



<a id="post-evaluation"></a>

#### Post-evaluation

The `dvf.inversion` function returns post-evaluation measures as optional
output. The measures are IC residual magnitude percentiles at each
iteration step.

Certain regions of the study domain are invalid for IC residual
calculations. These regions cover points that are either mapped outside the
domain by the forward DVF, or lie outside the deformed reference
domain. Function `dvf.maskDomain` calculates the valid domain.



<a id="gpu"></a>

#### GPU utilization

All `dvf.*` functions are GPU-enabled, supported by the MATLAB Parallel
Computing Toolbox. When a GPU is available, GPU computation is invoked by
casting the input data arrays into `gpuArray` objects.  For example:

    V = dvf.inversion(gpuArray(U)) 
	
The output data arrays are returned as `gpuArray` objects.



<a id="system-reqs"></a>

### System environment

The repository code was developed and tested on MATLAB R2018b. The
following functions in MATLAB toolboxes are used:

-   `imresize`, `imresize3`, `imdilate`, `imerode`, and `imclose` (Image
    Processing Toolbox, version 10.3);
-   `prctile` (Statistics and Machine Learning Toolbox, version 11.4); and
-   `gpuArray` (Parallel Computing Toolbox, version 6.13).

GPU computation was tested on an NVIDIA Quantum TXR113-1000R GPU and CUDA
10.0 drivers.



<a id="contributors"></a>

### Contributors

-   *Design, development, and testing:*  
    Alexandros-Stavros Iliopoulos, Abhishek Dubey, and Xiaobai Sun  
    Department of Computer Science, Duke University

-   *Consulting on medical CT image analysis and applications:*  
    Lei Ren and Fang-Fang Yin  
    Department of Radiation Oncology, Duke University Medical School
