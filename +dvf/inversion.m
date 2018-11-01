function [G, Mu, MaskCtrl, prctRGIter, prctRFIter] = inversion (F, varargin)
%
% INVERSION - Iterative inversion of 2D/3D deformation vector field with
%             bi-residual feedback control
%
% ----------------------------------------------------------------------
%
% Copyright (C) 2018, Department of Computer Science, Duke University
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% ----------------------------------------------------------------------
%   
% SYNTAX
%
%   G = INVERSION( F )
%   G = INVERSION( F, 'Name', Value, ... )
%   [G, MU, MASKCTRL] = INVERSION( F, ... )
%   [G, MU, MASKCTRL, PRCT_RG_K, PRCT_RF_K] = INVERSION( F, ... )
%
% INPUT
%
%   F           Forward displacement field      [NX x NY x NZ x 3]
%               (in voxel-length units)
%
% OUTPUT
%
%   G           Inverse displacement field      [NX x NY x NZ x 3]
%               estimate
%   MU          Feedback control parameter      [K-vector | NX x NY x NZ]
%               value(s)
%   MASKCTRL    Controllable region mask        [NX x NY x NZ] (logical)
%               (false over points where no
%               feasible control parameter value
%               exists for convergence)
%   PRCT_RG_K   Step-wise percentile measures   [5 x (K+1)]
%               of study IC residual
%   PRCT_RF_K   Step-wise percentile measures   [5 x (K+1)]
%               of reference IC residual
%
%       * All scalar/vector field array sizes in the documentation refer to
%         the case of 3D DVFs.  In the case of 2D DVFs, the sizes are
%         amended accordingly: 
%         - [NX x NY x NZ x 3] --> [NX x NY x 2]
%         - [NX x NY x NZ]     --> [NX x NY]
%         - etc...
%       * K is the total number of iteration steps
%       * PRCT_RG_K and PRCT_RF_K contain the 50th, 90th, 95th, 98th,
%         and 100th percentiles of IC residual magnitudes per step
%
% OPTIONS
%
%   'Control'           [string | A-vector | NXxNYxNZ]  {'midrange'}
%
%       Feedback control parameter scheme or value(s).  Numeric input
%       indicates fixed constant (scalar), alternating (A-vector), or
%       spatially-variant (NX x NY x NZ) control values.  Otherwise, one of
%       the following string IDs may be input to specify an adaptive
%       control value scheme [1]:
%
%       * 'midrange': local midrange values;
%       * 'optimal': pointwise optimal values (sensitive to initial guess);
%       * 'alternating': alternating values (midrange percentiles); or
%       * 'global-midrange': global midrange value.
%
%       If spatially variant control is used, MU_k(x) = MU(x + G_k(x)),
%       where MU is the control parameter field.
%
%   'Scale'             [scalar | 2-vector]             {[0.5 1.0]}
%
%       Domain rescaling factor(s) for multi-scale iteration.  A value of s
%       means scaling each dimension by s.
%
%   'NumIteration'      [scalar | 2-vector]             {10}
%
%       Max number of iteration steps per scale.
%
%   'StopThreshold'     [scalar | 2-vector]             {1e-6}
%
%       Error-magnitude threshold for iteration termination per scale.
%       (Magnitude measures are relative to voxel size in current scale.)
%
%   'AccThreshold'      [scalar]                        {1}
%
%       Error-bound threshold for acceleration stage initiation.  The
%       acceleration stage is triggered locally when the local error
%       magnitude bound is below the threshold.  Non-positive values
%       disable the acceleration stage.
%
%   'InitialValue'      [NX x NY x NZ x 3]              {zeros(NX,NY,NZ,3)}
%
%       Initial inverse DVF used in the iteration.
%
%   'Mask'              [NX x NY x NZ]                  {dvf.maskDomain(F)}
%
%       Mask over spatial domain of interest.  Percentile measures (IC
%       residuals, control values, etc) are defined only over the masked
%       region.
%
% DESCRIPTION
%
%   INVERSION iteratively computes an estimate of the inverse deformation G
%   given a forward deformation field F.  Inverse consistency (IC) residual
%   fields between F and G iterates are used as feedback in the iteration,
%   per the two-phase framework described and analyzed in Refs. [1,2].
%
%   The first phase of the iteration proceeds as [1]
%
%       G_k+1(x) = G_k(x) - (1 - MU_k(x)) * RG_k(x) ,
%
%   where
%
%       RG_k(x) = G_k(x) + F(x + G_k(x)) 
%
%   is the study IC residual at the k-th iteration step, and MU_k(x) is the
%   feedback control parameter value at x.
%
%   The IC residual RG is related to the unknown inversion error E by [1]
%
%       E_k(x) = JF(xi) * RG_k(x) ,
%
%   where JF is the Jacobian of the forward transformation, xi lies between
%   (x + G_k(x)) and (x + G*(x)), and G* is the (unknown) true inverse DVF.
%
%   The feedback control parameter MU may be a global constant, a set of
%   alternating global values throughout the iteration, or a scalar field
%   of spatially variant values.  By default, INVERSION uses spatially
%   variant control with adaptive local midrange values [1,2], to ensure
%   stable, linear-rate convergence towards the inverse.
%
%   Once the error has been made sufficiently small, ||E_k(x)||_oo <= 1,
%   the iteration proceeds to the second phase to accelerate convergence
%   [2]:
%
%       G_k+1(x) = G_k(x) - RF_k(x + G_k(x)) ,
%
%   where
%
%       RF_k(x) = F(x) + G_k(x + F(x))
%
%   is the reference IC residual.  Since the error is unknown, an upper
%   bound is used to guarantee the phase-2 condition.  In this phase, the
%   iteration converges quadratically.  This acceleration phase is referred
%   as an implicit Newton iteration, since it does not entail explicit
%   inversion of the Newton matrix (which suffers from several numerical
%   issues).
%
%   INVERSION uses a multi-scale iteration to facilitate phase integration
%   and transition.  By default, the iteration starts at 1/8-resolution
%   (1/4 in the 2D case) and then continues at full resolution.
%
% GPU SUPPORT
%
%   Inversion iteration computations are performed on the GPU if the
%   forward DVF (F) is a gpuArray object.
%
% NOTES
%
%   IC residual magnitude percentiles are only calculated if they are
%   specified by the user as output.  These calculations may substantially
%   increase the execution time per iteration step.
%
%   Several options that allow more fined-tuned control over the inversion
%   iteration computations have been removed from the INVERSION function
%   interface.  The interested user may easily re-enable them in the source
%   code.
%
% DEPENDENCIES
%
%   dvf.jacobian
%   dvf.eigJacobian
%   dvf.sizeVf
%   dvf.feedbackControlVal
%   dvf.coorDisplace
%   dvf.icResidual
%   dvf.interpnVf
%   dvf.rescale
%   dvf.maskDomain
%   util.ndgridGpu
%   util.parseOptArgs
%
% REFERENCES
%
%   [1] A. Dubey*, A.-S. Iliopoulos*, X. Sun, F.-F. Yin, and L. Ren,
%   "Iterative inversion of deformation vector fields with feedback
%   control," Medical Physics, vol. 45, no. 7, pp. 3147-3160, 2018.
%   DOI: 10.1002/mp.12962.
%
%   [2] A. Dubey, A.-S. Iliopoulos, X. Sun, F.-F. Yin, and L. Ren,
%   "Symmetric completion of deformation registration via bi-residual
%   inversion." In preparation.
%
%
% See also      dvf.feedbackControlVal, dvf.icResidual
%
    
    
    %% PARAMETERS
    
    % spatial domain size and dimensionality of DVF
    [szDom, dim] = dvf.sizeVf( F );
    if ~( ((dim == 3) && (length(szDom) == 3)) || ...
          ((dim == 2) && (length(szDom) == 2)) )
        error( [mfilename ':InvalidDimensionality'], ...
               ['Input DVF must be a NXxNYxNZx3 array (3D case)' ...
                ' or a NXxNYx2 array (2D case)'] );
    end
        
    % options and default values
    opt.control       = 'midrange';
    opt.scale         = [0.5, 1.0];
    opt.numIteration  = 10;
    opt.stopThreshold = 1e-6;
    opt.accThreshold  = 1;
    opt.InitialValue  = zeros( size(F), 'like', F );
    opt.Mask          = dvf.maskDomain( F );
    
    % static parameters
    accAsync          = true;              % asunchronous acceleration stage?
    interp            = 'linear';          % interpolation method
    prcIcMeasure      = [50 90 95 98 100]; % ICR percentile ranks
    icrFullResolution = true;              % calculate ICR in full resolution?
    windowControl     = 0.05 * min(szDom); % max local-control window size
    fallbackControl   = 0.8;               % mu value at incontrollable points
    prcMidAlternating = [98 50];           % alternating control percentiles
    prcErrorBound     = 99;                % Jacobian norm percentile rank
    valExtrap         = NaN;               % out-of-bounds sentinel value
    
    % parse optional arguments
    opt = util.parseOptArgs( opt, varargin{:} );
    
    % optional computation flags
    flag.prctR = (nargout > 3);
    
    % disable multi-scale iteration with more than 2 scales
    if length( opt.scale ) > 2
        error( [mfilename ':TooManyScales'], ...
               ['Multi-scale iteration with more than 2 scales is currently ' ...
                'disabled.'] );
    end
    
    
    %% INITIALIZATION
    
    % inverse estimate initialization
    G = opt.InitialValue;
    
    % if a domain mask was input as gpuArray, copy back to main DRAM (the
    % IMRESIZE function is quite problematic with gpuArrays)
    if isa( opt.Mask, 'gpuArray' )
        opt.Mask = gather( opt.Mask );
    end
    
    % expand termination criteria to match #scales
    if isscalar( opt.numIteration )
        opt.numIteration = repmat( opt.numIteration, [1, length(opt.scale)] );
    end
    if isscalar( opt.stopThreshold )
        opt.stopThreshold = repmat( opt.stopThreshold, [1, length(opt.scale)] );
    end
    
    % preallocate space for optional-output arrays
    prctRGIter = nan( length(prcIcMeasure), sum(opt.numIteration) + 1, ...
                      'like', F );
    prctRFIter = nan( size(prctRGIter), 'like', F );
    
    
    %% PRE-PROCESSING
    
    [Mu, bndNormJF, szDomScl, MaskCtrl] = ...
        preprocessing( F, opt.control, opt.scale, opt.Mask, windowControl, ...
                       fallbackControl, prcErrorBound, prcMidAlternating, ...
                       opt.numIteration, dim );
    
    
    %% DVF INVERSION ITERATION
    
    % "total" iteration counter (across all scales)
    ks = 0;
    
    for s = 1 : length(opt.scale)       % ======== each scale
        
        % rescale forward and (current) inverse DVFs, as well as the
        % domain mask and control parameter map
        [Fs, G, MaskDomainS, MaskCtrlS, MuS] = ...
            prepareScale( F, G, opt.Mask, MaskCtrl, Mu, szDomScl{s}, interp );
        
        % rescaled reference domain grid coordinates
        grid      = cell( 1, dim );
        [grid{:}] = util.ndgridGpu( szDomScl{s}, isa( F, 'gpuArray' ) );
        
        % perform iteration in current scale
        [G, k, prctRGIter(:,ks+(1:opt.numIteration(s))), ...
               prctRFIter(:,ks+(1:opt.numIteration(s)))] = ...
            iterationCore( Fs, G, MuS, MaskDomainS, grid, MaskCtrlS, ...
                           opt.stopThreshold(s), opt.accThreshold, ...
                           bndNormJF, opt.numIteration(s), interp, ...
                           valExtrap, opt.scale(s)*windowControl, ...
                           fallbackControl, prcIcMeasure, szDom, accAsync, ...
                           flag.prctR, icrFullResolution );
        
        % scale back to original domain size
        G = dvf.rescale( G, szDom );
        
        % update "total" iteration counter
        ks = ks + k;
        
    end  % for (s)
    
    
    %% TERMINATION
    
    % total number of iterations (+1 to account for initial error/residual)
    ks = ks + 1;
    
    % final IC residual measures
    if flag.prctR
        [~, ~, ~, ~, prctRGIter(:,ks), prctRFIter(:,ks)] = ...
            icResidualCalc( F, G, interp, valExtrap, {}, opt.Mask, ...
                            prcIcMeasure, flag.prctR );
    end  % if (final IC residual measures)
    

end



%% ==================================================
%  PRE-PROCESSING

function [Mu, bndJFNorm, szDomScl, MaskCtrl, JF, LambdaF] = ...
        preprocessing (F, control, scales, Mask, wndCtrl, valCtrlFallback, ...
                       prcErrBnd, prcMidAlt, numIter, dim)
% IN    F               forward DVF             [NX x NY x NZ x 3]
%       control         feedback control param. [<numeric> | string]
%       tauAcc          acceleration threshold  [scalar]
%       scales          iteration resol. scales [S-vector]
%       Mask            domain of interest mask [NX x NY] (logical)
%       wndCtrl         feedback control n/hood [scalar | WX x WY x WZ]
%       valCtrlFallback control fallback value  [scalar]
%       prcErrBnd       % for Jacob. norm bound [scalar]
%       prcMidAlt       alternating control %s  [A-vector]
%       dim             domain dimensionality   [2 | 3]
%       numIter         # of iteration steps    [S-vector]
% OUT   Mu              control parameter values[K-vector | NX x NY x NZ]
%       bndJFNorm       ||JF(x)||_oo bound      [scalar]
%       szDomScl        scaled domain sizes     [S-cell(3-vector)]
%       MaskCtrl        controllable domain mask[NX x NY x NZ] (logical)
%       JF              deformation Jacobian    [NX x NY x NZ x 3 x 3]
%       LambdaF         Jacobian eigenvalues    [NX x NY x NZ x 3]
    
    % forward deformation Jacobian
    JF = dvf.jacobian( F );
    
    % non-invertible subdomain
    MaskNonInv = noninvertibleDomain( JF );
    
    % feedback control values set-up
    if ischar( control )                % ---- NTDC-adaptive control
        
        % Jacobian spectrum calculation
        LambdaF = dvf.eigJacobian( JF );
        
        % local control neighborhood (if adaptive control neighborhoods are used,
        % set the trivial neighborhood---local calculations will be done
        % during the iteration rather than preprocessing)
        if isscalar( wndCtrl )
            wndCtrl = 1;
        end
        
        % feedback control parameter values set-up
        [Mu, MaskCtrl] = dvf.feedbackControlVal( ...
            LambdaF, control, 'Mask', Mask, ...
            'WindowMidrange', wndCtrl, ...
            'FallbackValue', valCtrlFallback, ...
            'PrcMidAlternating', prcMidAlt );
        
        % controllable region is only relevant within the user-specified domain
        MaskCtrl = MaskCtrl & Mask;
        
    else                                % ---- user-input control
        
        Mu       = control;
        MaskCtrl = Mask;
        
    end  % if (adaptive/prefixed control?)
    
    % non-invertible & incontrollable domain warning
    if any( MaskNonInv, 'all' )
        warning( [mfilename ':NonInvertibleRegion'], ...
                 'Input DVF is non-invertible at %d/%d points', ...
                 nnz(MaskNonInv), numel(MaskNonInv) );
    end
    if any( ~MaskCtrl, 'all' )
        warning( [mfilename ':IncontrollableRegion'], ...
                 'The inversion error is incontrollable at %d/%d points', ...
                 nnz(~MaskCtrl), numel(MaskCtrl) );
    end
    MaskCtrl = MaskCtrl & ~MaskNonInv;
    
    % in case of constant/alternating control values, replicate for convenience
    if isvector( Mu )
        Mu = repmat( Mu(:), [ceil( max(numIter) * length(Mu) ), 1] );
    end
    
    % "soft" upper bound of Jacobian oo-norm over controllable/masked domain
    bndJFNorm = max( sum( abs(JF), dim+2 ), [], dim+1 );
    bndJFNorm = prctile( bndJFNorm(MaskCtrl), prcErrBnd );
    
    % turn scale factors to scaled-domain sizes to avoid floor/ceiling issues
    % when scaling down and then up with IMRESIZE (e.g. 11->6->12)
    [szDom, ~] = dvf.sizeVf( F );
    szDomScl   = arrayfun( @(s) ceil( szDom * s ), scales, ...
                           'UniformOutput', false );
    
end



%% ==================================================
%  INVERSION ITERATION CORE

function [G, k, prctRGIter, prctRFIter] = ...
        iterationCore (F, G, Mu, MaskDom, grid, MaskCtrl, tauStop, tauAcc, ...
                       bndNormJF, numIter, interp, extrap, radiusMuMax, ...
                       muValOob, prcPop, szDomFull, flagAccAsync, flagPrctR, ...
                       flagIcFullRes)
% IN    F               forward DVF             [NX x NY x NZ x 3]
%       G               inverse DVF estimate    [NX x NY x NZ x 3]
%       Mu              feedback control values [K-vector | NX x NY x NZ]
%       MaskDom         spatial domain mask     [NX x NY x NZ] (logical)
%       MaskCtrl        controllable region mask[NX x NY x NZ] (logical)
%       grid            spatial grid coordinates[3-cell(NX x NY x NZ)]
%       tauStop         termination threshold   [scalar]
%       tauAcc          acceleration threshold  [scalar]
%       bndNormJF       Jacobian oo-norm bound  [scalar]
%       numIter         # of iteration steps    [scalar]
%       interp          interpolation method    [string]
%       extrap          extrapolation value     [scalar]
%       radiusMuMax     max mu window radius    [scalar]
%       muValOob        out-of-bounds mu value  [scalar]
%       prcPop          ICR percentile ranks    [P-vector]
%       szDomFull       domain size in full res.[3-vector]
%       flagAccAsync    asynchronous accel.?    [boolean]
%       flagPrctR       calculate ICR prctiles? [boolean]
%       flagIcFullRes   full-res. ICR prctiles? [boolean]
% OUT   G               next inverse DVF estim. [NX x NY x NZ x 3]
%       k               # of taken iter. steps  [scalar]
%       prctRGIter      stepwise RG percentiles [P x K]
%       prctRFIter      stepwise RF percentiles [P x K]

    % ========== INITIALIZATION
    
    % preallocate space for optional output arrays
    if flagPrctR
        prctRGIter = nan( length(prcPop), numIter );
        prctRFIter = nan( length(prcPop), numIter );
    end
    
    % initialize full-resolution data if required for IC residual
    % magnitude percentile calculations
    if flagIcFullRes && ~isequal( dvf.sizeVf(G), szDomFull )
        MaskDomS = resizeMask( MaskDom, szDomFull );
        FS       = dvf.rescale( F, szDomFull );
    end
    
    % ========== ITERATION
    
    for k = 1 : numIter
        
        % deformed reference grid coordinates
        gridG      = cell( size(grid) );
        [gridG{:}] = dvf.coorDisplace( G, grid{:} );
        
        % IC residual vector fields, magnitudes, and percentiles
        [RG, RF, RGNorm, RFNorm, prctRGIter(:,k), prctRFIter(:,k)] = ...
            icResidualCalc( F, G, interp, extrap, gridG, MaskDom, ...
                            prcPop, flagPrctR );
        
        % IC residual percentiles calculation in full resolution?
        if flagIcFullRes && ~isequal( dvf.sizeVf(G), szDomFull )
            [~, ~, ~, ~, prctRGIter(:,k), prctRFIter(:,k)] = ...
                icResidualCalc( FS, dvf.rescale( G, szDomFull ), ...
                                interp, extrap, {}, MaskDomS, ...
                                prcPop, flagPrctR );
        end  % if (full-resolution ICR percentiles)
        
        % ||E(x)||_oo  = ||JF(xi) * RG(x)||_oo
        %             <= ||JF(xi)||_oo * ||RG(x)||_oo
        %             <= max_z||JF(z)||_oo * ||RG(x)||_oo
        %
        % -->  ||RG(x)||_oo >= ||E(x)||_oo / max_z||JF(x)||_oo
        
        % termination criterion check
        if all( RGNorm(MaskCtrl) < tauStop/bndNormJF, 'all' ) ...
                || (max(RFNorm(MaskCtrl)) < tauStop)
            break
        end
        
        % local acceleration domain mask
        MaskAcc = (RGNorm < tauAcc / bndNormJF);
        if ~flagAccAsync                % synchronous (global) acceleration?
            MaskAcc = all( MaskAcc(Mask), 'all' );
        end
        
        % neighborhood radius and mask of points with small enough radius of
        % neighborhood for adaptive local control---the rest of the points
        % are assigned a "global" control value (over the sub-domain that
        % covers them)
        if isscalar( radiusMuMax )
            radiusMu = min( [max( RFNorm(MaskCtrl) ), ...
                             max( RGNorm(MaskCtrl) ) * bndNormJF, ...
                             radiusMuMax] );
        else
            radiusMu  = 0;
        end
        
        % update inverse DVF iterate
        G = G - computeFeedback( RG, RF, Mu, MaskAcc, k, interp, extrap, ...
                                 gridG, radiusMu, muValOob );
        
    end  % for (k; iteration steps)
    
end



%% ==================================================
%  FEEDBACK FIELD COMPUTATION

function Fb = computeFeedback (RG, RF, Mu, MaskAcc, k, interp, ...
                               extrap, gridG, radiusMu, muValOob)
% IN    RG              study IC residual       [NX x NY x NZ x 3]
%       RF              reference IC residual   [NX x NY x NZ x 3]
%       Mu              control parameter       [scalar |vector | NX x NY x NZ]
%       MaskAcc         acceleration domain mask[NX x NY x NZ] (logical)
%       k               iteration step #        [scalar]
%       interp          interpolation method    [string]
%       extrap          extrapolation value     [scalar]
%       gridG           G-deformed domain grid  [3-cell(NX x NY x NZ)]
%       radiusMu        max-error radius        [scalar]
%       muValOob        out-of-bounds mu value  [scalar]
% OUT   Fb              feedback field iterate  [NX x NY x NZ x 3]

    % domain dimensionality
    [~, dim] = dvf.sizeVf( RG );
    
    % ========== SCALAR-CONTROL FEEDBACK
    
    % skip scalar-control calculations if the whole domain is in the
    % acceleration phase
    if ~all( MaskAcc, 'all' )           % -------- some scalar control
    
        if isvector( Mu )               % ---- uniform value
            
            Mu = Mu(k);
            
        else                            % ---- spatially variant values
            
            % look for most conservative control value in adaptively-sized
            % neighborhoods? (IMDILATE only supports UINT8 and LOGICAL input on
            % the GPU...)
            if radiusMu > 0
                windowMu = true( repmat( ceil(2*radiusMu), [1 dim] ) );
                Mu       = imdilate( gather(Mu), windowMu, 'same' );
                if isa( RG, 'gpuArray' )
                    Mu = gpuArray( Mu );
                end
            end
            
            % control values map based on current-iterate deformed grid
            % MU_k(x) = MU(x + G_k(x))
            Mu = interpn( Mu, gridG{:}, interp, extrap );
            Mu( isnan(Mu) ) = muValOob;
            
        end  % if (global/variant control)
        
        % scalar-control feedback field
        FbSC = (1 - Mu) .* RG;
        
    else                                % -------- no more scalar control
        
        FbSC = 0;
        
    end  % if (don't skip straight to implicit Newton step)
    
    % ========== IMPLICIT NEWTON FEEDBACK
    
    % RF(y), y = x + G(x) is used directly as feedback
    
    % ========== FEEDBACK FIELD MASKING
    
    Fb = (FbSC .* ~MaskAcc) + (RF .* MaskAcc);
    
    
end



%% ==================================================
%  IC RESIDUAL CALCULATIONS

function [RG, RF, RGNorm, RFNorm, prctRG, prctRF] = ...
        icResidualCalc (F, G, interp, extrap, gridG, Mask, prcPop, flagPrct)
% IN    F               forward DVF             [NX x NY x NZ x 3]
%       G               inverse DVF estimate    [NX x NY x NZ x 3]
%       interp          interpolation method    [string]
%       extrap          extrapolation value     [scalar]
%       gridG           G-deformed grid coors.  [3-cell(NX x NY x NZ)]
%       Mask            domain mask for prctiles[NX x NY x NZ] (logical)
%       prcPop          population measure %s   [P-vector]
%       flagPrct        calculate ICR prctiles? [boolean]
% OUT   RG              study IC residual       [NX x NY x NZ x 3]
%       RF              reference IC residual   [NX x NY x NZ x 3]
%       RGNorm          RG infinity norms       [NX x NY x NZ]
%       RFNorm          RF infinity norms       [NX x NY x NZ]
%       prctRG          RG percentiles          [P-vector]
%       prctRF          RF percentiles          [P-vector]

    % domain dimensionality
    [~, dim] = dvf.sizeVf( F );
    
    % deformed-grid coordinates
    if isempty( gridG )
        gridG      = cell( 1, dim );
        [gridG{:}] = dvf.coorDisplace( G );
    end
    
    % IC residual:  RG(x) = G(x) + F(x + G(x))
    RG = dvf.icResidual( G, F, interp, extrap, gridG{:} );
    
    % IC residual:  RF(y) = F(y) + G(y + F(y))          [y = x + G(x)]
    %                     = (F(x+G(x)) + G(x)) + G(x+G(x)+F(x+G(x))) - G(x)
    %                     = RG(x) + G(x + RG(x)) - G(x)
    % *** RF(x + G(x)) can be numerically unstable (in part because it
    %     involves non-small deformations when RF is used as feedback in
    %     the implicit Newton iteration), plus it is inefficient compared
    %     to the formulation below (involving an additional interpolation
    %     step)
    gridRG      = cell( size(gridG) );
    [gridRG{:}] = dvf.coorDisplace( RG );
    RF          = RG + (dvf.interpnVf( G, gridRG, interp, extrap ) - G);
    
    % IC residual magnitude percentile measures
    prctRG = icrMgnPercentiles( RG, Mask, prcPop, flagPrct );
    prctRF = icrMgnPercentiles( RF, Mask, prcPop, flagPrct );
    
    % IC residual field infinity-norms
    RGNorm = max( abs(RG), [], dim+1 );
    RFNorm = max( abs(RF), [], dim+1 );
    
    % deal with potential NaN/Inf values (out-of-bounds interpolation)
    [RG, RGNorm] = sanitizeIcResidual( RG, RGNorm );
    [RF, RFNorm] = sanitizeIcResidual( RF, RFNorm );
    
end



%% ==================================================
%  FIELD & MASK RESCALING PRIOR TO ITERATION CORE

function [FS, GS, MaskDomainS, MaskCtrlS, MuS] = ...
        prepareScale (F, G, MaskDomain, MaskCtrl, Mu, szDomS, interp)
% IN    F               forward DVF             [NX x NY x NZ x 3]
%       G               inverse DVF estimate    [NX x NY x NZ x 3]
%       MaskDomain      valid domain mask       [NX x NY x NZ] (logical)
%       MaskCtrl        controllable region mask[NX x NY x NZ] (logical)
%       Mu              feedback control values [K-vector | NX x NY x NZ]
%       szDomS          rescaled domain size    [SX SY SZ]
%       interp          DVF interpolation method[string]
% OUT   FS              rescaled F              [SX x SY x SZ x 3]
%       GS              rescaled G              [SX x SY x SZ x 3]
%       MaskDomainS     rescaled MaskDomain     [SX x SY x SZ] (logical)
%       MaskCtrlS       rescaled MaskCtrlS      [SX x SY x SZ] (logical)
%       MuS             rescaled MuS (if needed)[K-vector | SX x SY x SZ]
    
    % DVFs
    FS = dvf.rescale( F, szDomS, interp );
    GS = dvf.rescale( G, szDomS, interp );
    
    % masks
    MaskDomainS = resizeMask( MaskDomain, szDomS );
    MaskCtrlS   = resizeMask( MaskCtrl  , szDomS );
    
    % feedback control values
    MuS = resizeMu( Mu, szDomS );
    
end



%% ==================================================
%  DOMAIN MASK RESIZING

function MaskS = resizeMask (Mask, szDomS)
% IN    Mask            domain mask             [NX x NY x NZ] (logical)
%       szDomS          resized domain size     [SX SY SZ]
% OUT   MaskS           resized domain mask     [SX x SY x SZ] (logical)
    
    % domain dimensionality
    dim = length( szDomS );
    
    % resize domain mask
    switch dim
        
      case 2                            % ---- 2D
        MaskS = logical( imresize( Mask, szDomS, 'bilinear' ) );
        
      case 3                            % ---- 3D
        MaskS = logical( imresize3( uint8(Mask), szDomS, 'linear' ) );
        % (IMRESIZE3 does not work with logical data)
        
    end  % switch (2D/3D case)
    
end



%% ==================================================
%  FEEDBACK CONTROL VALUE MAP RESIZING

function MuS = resizeMu (Mu, szDomS)
% IN    Mu              feedback control values [NX x NY x NZ | K-vector]
%       szDomS          resized domain size     [SX SY SZ]
% OUT   MuS             resized Mu values       [SX x SY x SZ]
    
    % if spatially uniform values are used, do nothing
    if isvector( Mu )
        MuS = Mu;
        return
    end
    
    % domain dimensionality
    dim = length( szDomS );
    
    % resize feedback control parameter value map
    switch dim
        
      case 2                            % ---- 2D
        MuS = imresize( Mu, szDomS, 'bilinear' );
        
      case 3                            % ---- 3D
        MuS = imresize3( Mu, szDomS, 'linear' );
        
    end  % switch (2D/3D case)
    
end



%% ==================================================
%  SANITIZE NAN/INF INTERPOLATION VALUES IN IC RESIDUAL

function [R, RMgn] = sanitizeIcResidual (R, RMgn)
%       R               IC residual vectorfield [NX x NY x NZ x 3]
%       RMgn            IC residual magnitudes  [NX x NY x NZ]
    
    % domain dimensionality
    dim = ndims( RMgn );
    
    % spatial mask of "bad" interpolation values (NaN/Inf)
    MaskBad = isnan( RMgn ) | isinf( RMgn );
    
    % set invalid IC residual magnitudes equal to max IC residual (in
    % order to not skew statistics)
    RMgn(MaskBad) = max( RMgn(~MaskBad) );
    
    % set invalid IC residual vectors equal to 0 (no feedback control)
    R = reshape( R, [numel(MaskBad), dim] );
    R(MaskBad,:) = 0;
    R = reshape( R, [size(MaskBad), dim] );
    
end



%% ==================================================
%  IC RESIDUAL MAGNITUDE PERCENTILES

function prctR = icrMgnPercentiles (R, Mask, prcPop, flag)
% IN    R               IC residual             [NX x NY x NZ x 3]
%       Mask            spatial domain mask     [NX x NY x NZ] (logical)
%       prcPop          ICR percentile ranks    [P-vector]
%       flag            compute ICR percentiles?[boolean]
% OUT   prctR           ICR magnitude prctiles  [P-vector]
    
    if flag                             % ---- compute ICR percentiles
        dim       = ndims( Mask );
        RMgn      = sqrt( sum( R.^2, dim+1 ) );
        [~, RMgn] = sanitizeIcResidual( R, RMgn );
        prctR     = prctile( RMgn(Mask), prcPop, 'all' );
    else                                % ---- do not compute anything
        prctR = NaN;
    end  % if
    
end



%% ==================================================
%  NON-INVERTIBLE DOMAIN MASK

function MaskNonInv = noninvertibleDomain (JF)
% IN    JF              deformation Jacobian    [NX x NY x NZ x 3 x 3]
% OUT   MaskNonInv      noninvertible dom. mask [NX x NY x NZ] (logical)
    
    % numerical zero threshold
    tau = 1e-5;
    
    % Jacobian determinant (explicit calculation)
    % *** WARNING numerically sensitive
    DetJF = JF(:,:,:,1,1) * (JF(:,:,:,2,2) * JF(:,:,:,3,3) - ...
                             JF(:,:,:,2,3) * JF(:,:,:,3,2)) - ...
            JF(:,:,:,1,2) * (JF(:,:,:,2,1) * JF(:,:,:,3,3) - ...
                             JF(:,:,:,2,3) * JF(:,:,:,3,1)) + ...
            JF(:,:,:,1,3) * (JF(:,:,:,2,1) * JF(:,:,:,3,2) - ...
                             JF(:,:,:,2,2) * JF(:,:,:,3,1));
    
    % DVF is non-invertible at all points where the Jacobian is singular
    MaskNonInv = abs(DetJF) <= tau;
    
end



%%------------------------------------------------------------
%
% AUTHORS
%
%   Alexandros-Stavros Iliopoulos       ailiop@cs.duke.edu
%
% VERSION
%
%   1.0.0 - October 31, 2018
%
% ------------------------------------------------------------