# idvf

Iterative Inversion of Deformation Vector Field with Adaptive Bi-residual
Feedback Control



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
denoted by <img alt="$u$"
src="https://rawgit.com/ailiop/idvf/None/svgs/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode"
align="middle" width="9.41028pt" height="14.15535pt"/>. The DVF describes a
(non-linear) mapping from a reference image <img alt="$I_r$"
src="https://rawgit.com/ailiop/idvf/None/svgs/68a1e29bb03da1906d9220ebd533ad09.svg?invert_in_darkmode"
align="middle" width="13.683615pt" height="22.46574pt"/> onto a study image
<img alt="$I_s$"
src="https://rawgit.com/ailiop/idvf/None/svgs/0096b70ae460e0de5f398f022769a8f7.svg?invert_in_darkmode"
align="middle" width="13.430505pt" height="22.46574pt"/> via point-wise
displacement, i.e., <p align="center"><img alt="$$ &#10;I_r(y) = I_s(y +
u(y)), &#10;$$"
src="https://rawgit.com/ailiop/idvf/None/svgs/1aaeda0abb91834a3b58c6321b20a8af.svg?invert_in_darkmode"
align="middle" width="149.046975pt" height="16.438356pt"/></p> at all
locations <img alt="$y$"
src="https://rawgit.com/ailiop/idvf/None/svgs/deceeaf6940a8c7a5a02373728002b0f.svg?invert_in_darkmode"
align="middle" width="8.6493pt" height="14.15535pt"/> in the reference
domain. We aim at computing the inverse DVF <img alt="$v$"
src="https://rawgit.com/ailiop/idvf/None/svgs/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode"
align="middle" width="8.55789pt" height="14.15535pt"/> such that <p
align="center"><img alt="$$ &#10;I_s(x) = I_r(x + v(x)), &#10;$$"
src="https://rawgit.com/ailiop/idvf/None/svgs/e77b118cb0d39454a3696351c4aac28e.svg?invert_in_darkmode"
align="middle" width="150.43182pt" height="16.438356pt"/></p> at all <img
alt="$x$"
src="https://rawgit.com/ailiop/idvf/None/svgs/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode"
align="middle" width="9.3951pt" height="14.15535pt"/> in the study domain.


<a id="ic-residuals"></a>

#### Two inverse consistency residuals

If two DVFs <img alt="$u$"
src="https://rawgit.com/ailiop/idvf/None/svgs/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode"
align="middle" width="9.41028pt" height="14.15535pt"/> and <img alt="$v$"
src="https://rawgit.com/ailiop/idvf/None/svgs/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode"
align="middle" width="8.55789pt" height="14.15535pt"/> are inverse of each
other, they must meet the *inverse consistency* (IC) condition. We define
the study IC residual, <p align="center"><img alt="$$&#10;s(x) = v(x) +
u(x + v(x)),&#10;$$"
src="https://rawgit.com/ailiop/idvf/None/svgs/6071c793b7bcaacacfa7b12b6173f200.svg?invert_in_darkmode"
align="middle" width="189.61965pt" height="16.438356pt"/></p> and the
reference IC residual, <p align="center"><img alt="$$&#10;r(y) = u(y) +
v(y + u(y)).&#10;$$"
src="https://rawgit.com/ailiop/idvf/None/svgs/c2f4b6419bb41288ad987829a8ccf10f.svg?invert_in_darkmode"
align="middle" width="187.65615pt" height="16.438356pt"/></p> The DVFs <img
alt="$u$"
src="https://rawgit.com/ailiop/idvf/None/svgs/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode"
align="middle" width="9.41028pt" height="14.15535pt"/> and <img alt="$v$"
src="https://rawgit.com/ailiop/idvf/None/svgs/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode"
align="middle" width="8.55789pt" height="14.15535pt"/> are
inverse-consistent if both residuals are zero. The residuals are related to
the inversion error, <img alt="$e(x) = v(x) -&#10;v_*(x)$"
src="https://rawgit.com/ailiop/idvf/None/svgs/768f11fd69defab2835e3445913d2f7a.svg?invert_in_darkmode"
align="middle" width="140.28729pt" height="24.6576pt"/>, where <img
alt="$v_*$"
src="https://rawgit.com/ailiop/idvf/None/svgs/62555e11fd1268ce81658e8d04041225.svg?invert_in_darkmode"
align="middle" width="14.703315pt" height="14.15535pt"/> is the true
inverse DVF. The IC residuals are computationally available, while the
inversion error is unknown. We use the IC residuals as feedback in
iterative DVF inversion, and exercise adaptive control over the feedback
for global convergence and local acceleration.


<a id="inversion-iteration"></a>

#### DVF inversion iteration with bi-residual feedback control

Due to the non-linear nature of the DVF mapping, one resorts to iterative
procedures for DVF inversion. Our procedure has two stages: preprocessing
and adaptive inversion iteration.

In the preprocessing stage, we compute the DVF Jacobians and their
eigenvalues over the spatial domain. We use this spectral information to
determine data-adaptive parameters for feedback control.

The inversion iteration proceeds in two phases. During phase one, at step
<img alt="$k+1$"
src="https://rawgit.com/ailiop/idvf/None/svgs/33359de825e43daa97171e27f6558ae9.svg?invert_in_darkmode"
align="middle" width="37.385865pt" height="22.83138pt"/>, we update the
current inverse estimate <img alt="$v_k$"
src="https://rawgit.com/ailiop/idvf/None/svgs/eaf0887cdc4cb5f8e69a7796f143c3eb.svg?invert_in_darkmode"
align="middle" width="15.23412pt" height="14.15535pt"/> using the study IC
residual <img alt="$s_k$"
src="https://rawgit.com/ailiop/idvf/None/svgs/59efeb0f4f5d484a9b8a404d5bdac544.svg?invert_in_darkmode"
align="middle" width="14.971605pt" height="14.15535pt"/> as feedback: <p
align="center"><img alt="$$&#10;v_{k+1}(x) = v_k(x) - (1 - \mu_k(x)) s_k(x)
,&#10;$$"
src="https://rawgit.com/ailiop/idvf/None/svgs/8a19510dca81766687d4e19803899701.svg?invert_in_darkmode"
align="middle" width="258.9345pt" height="16.438356pt"/></p> where <img
alt="$\mu$"
src="https://rawgit.com/ailiop/idvf/None/svgs/07617f9d8fe48b4a7b3f523d6730eef0.svg?invert_in_darkmode"
align="middle" width="9.90495pt" height="14.15535pt"/> is the feedback
control parameter. With the adaptive control schemes provided in this
repository, the iteration is guaranteed to converge globally, under certain
mild conditions [<a href="#medphys2018">1</a>].

Inversion errors can be estimated using the IC residuals and the DVF
spectral information. Once the errors are made sufficiently small, the
iteration is switched to phase two: <p align="center"><img
alt="$$&#10;v_{k+1}(x) = v_k(x) - r_k(x + v_k(x)),&#10;$$"
src="https://rawgit.com/ailiop/idvf/None/svgs/cfa2d7ff960254eb3b4e70cf32b5aec6.svg?invert_in_darkmode"
align="middle" width="235.70415pt" height="16.438356pt"/></p> where the
reference IC residual <img alt="$r_k$"
src="https://rawgit.com/ailiop/idvf/None/svgs/eed77c5296d3cd11c33cd86d1e14efef.svg?invert_in_darkmode"
align="middle" width="14.68236pt" height="14.15535pt"/> at displaced
locations is used as feedback. The local convergence rate is quadratic
during phase two. We may refer to this iteration as implicit Newton
iteration. However, phase-two iteration steps do not entail explicit
formation and inversion of the Newton matrix; the explicit Newton step
suffers from several numerical problems.

Phase transition and integration in the inversion iteration is enabled by a
multi-resolution scheme, elaborated in Reference [<a
href="#dukephd2018">2</a>].


<a id="references"></a>

#### References

<a id="medphys2018"></a>[1] A. Dubey*, A.S. Iliopoulos*, X. Sun, F.F. Yin,
and L. Ren, <a href="http://dx.doi.org/10.1002/mp.12962">Iterative
inversion of deformation vector fields with feedback control</a>,
*Medical Physics*, 45(7):3147-3160 (2018).

<a id="dukephd2018"></a>[2] A. Dubey, *Symmetric Completion of Deformable
Registration via Bi-residual Inversion*, PhD dissertation, Duke
University, Durham, NC, USA, 2018.


<a id="software"></a>

### Software description


<a id="inversion-function"></a>

#### DVF inversion

The main function is `dvf.inversion`. It takes a forward DVF <img alt="$u$"
src="https://rawgit.com/ailiop/idvf/None/svgs/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode"
align="middle" width="9.41028pt" height="14.15535pt"/> and returns the
inverse DVF <img alt="$v$"
src="https://rawgit.com/ailiop/idvf/None/svgs/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode"
align="middle" width="8.55789pt" height="14.15535pt"/>. A 3D DVF is
represented by a 4D array of size <img alt="$N_x \times N_y \times N_z
\times 3$"
src="https://rawgit.com/ailiop/idvf/None/svgs/e43ae8c45f6609ce1d55c352b19fb71b.svg?invert_in_darkmode"
align="middle" width="131.868165pt" height="22.46574pt"/>, in `single` or
`double` precision. Specifically, `U(i,j,k,1:3)` is the 3D forward
displacement vector <img alt="$u(x)$"
src="https://rawgit.com/ailiop/idvf/None/svgs/320b3450fd8b780975b68c70115439b3.svg?invert_in_darkmode"
align="middle" width="31.590735pt" height="24.6576pt"/> at voxel <img
alt="$x = (i,j,k)$"
src="https://rawgit.com/ailiop/idvf/None/svgs/80421b97532d6f6d46d108ac51b2848e.svg?invert_in_darkmode"
align="middle" width="80.245605pt" height="24.6576pt"/>. Two-dimensional
DVFs are represented similarly.

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
parameter <img alt="$\mu$"
src="https://rawgit.com/ailiop/idvf/None/svgs/07617f9d8fe48b4a7b3f523d6730eef0.svg?invert_in_darkmode"
align="middle" width="9.90495pt" height="14.15535pt"/> per the chosen
scheme.

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
casting the input data arrays into `gpuArray` objects; e.g., `V =
dvf.inversion(gpuArray(U))`. The output data arrays are returned as
`gpuArray` objects.


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
