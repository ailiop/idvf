function [Mu, MaskCtrl, P] = feedbackControlVal (Lambda, control, varargin)
%
% FEEDBACKCONTROLVAL - Feedback control values and error propagation
%                      factors for DVF inversion iteration 
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
%   [MU, MASKCTRL] = FEEDBACKCONTROLVAL( LAMBDA )
%   [MU, MASKCTRL] = FEEDBACKCONTROLVAL( LAMBDA, CONTROL )
%   [MU, MASKCTRL] = FEEDBACKCONTROLVAL( LAMBDA, CONTROL, 'Name', Value, ... )
%   [MU, MASKCTRL, P] = FEEDBACKCONTROLVAL( LAMBDA, CONTROL, ... )
%
% INPUT
%
%   LAMBDA      Deformation Jacobian eigenvalues[NX x NY (x NZ) x (2|3)]
%   CONTROL     Feedback control scheme         [string | vector | NX x NY (x NZ)]
%               * 'midrange'            (local midrange values)
%               * 'optimal'             (pointwise optimal values; sensitive)
%               * 'alternating'         (population-midrange alternating values)
%               * 'global-midrange'     (global midrange value)
%               * [scalar]              (constant global value)
%               * [M-vector]            (M alternating global values)
%               * [NX x NY x NZ]        (spatially variant values)
%               {default: 'midrange'}
%
% OUTPUT
%
%   MU          Feedback control values         [NX x NY (x NZ) | vector]
%   MASKCTRL    Controllable region mask        [NX x NY (x NZ)]
%   P           Local error propagation factors [NX x NY (x NZ)]
%
% OPTIONS
%
%   'Mask'              [NX x NY (x NZ)]                {true(NX,NY(,NZ))}
%
%       Logical array.  Mask over the spatial domain, indicating the domain
%       for global calculations ('alternating' and 'global-midrange'
%       control schemes).
%
%   'WindowMidrange'    [scalar | WX x WY (x WZ)]       {1}
%
%       Local neighborhood window size or mask for local midrange
%       control value search.
%       * [scalar]              square/cubic window side-length
%       * [WX x WY (x WZ)]      2D/3D neighborhood mask
%
%   'PrcMidAlternating' [A-vector]                      {[98 50]}
%
%       Percentile ranks for alternating feedback control.  The alternating
%       control values are set equal to the corresponding midrange value
%       percentiles.
%
%   'FallbackValue'     [scalar]                        {0.8}
%
%       Control parameter value to use at all points which lie outside the
%       input mask or at which the DVF controllability condition [1] is
%       violated.
%
% DESCRIPTION
%
%   FEEDBACKCONTROLVAL sets feedback control parameter values for the DVF
%   inversion iteration [1] based on the input control scheme, and
%   calculates the corresponding point-wise (infinitesimal) error
%   propagation factors.
%
%   If the eigenvalues array (LAMBDA) contains any NaN values, the
%   points/indices where NaNs are encountered are treated as "masked off"
%   of the spatial domain.  The feedback control parameter value (MU) at
%   these locations is set to the fallback value, and no error propagation
%   factors are calculated (NaN is returned for those locations in P).
%
% ASSUMPTION
%
%   The input eigenvalues refer to a real-valued Jacobian.  Therefore, in
%   the 3D case, the eigenvalues are either all real, or 1 is real and the
%   other 2 are a conjugate pair; in the 2D case, they are either both real
%   or a conjugate pair.
%
% SCHEMES
%
%   * 'midrange'
%
%       Local midrange control values [1,2].  The local neighborhood
%       size is defined by the WindowMidrange option.
%
%       MU(x) = 1 - gamma(N(x))
%
%       where gamma(x) = min_j{ Re(L_j(x)) / |L_j(x)|^2 }, N(x) is a
%       neighborhood around x, and L refers to the deformation Jacobian
%       eigenvalues.
%
%   * 'optimal'
%
%       Locally optimal control values [1].
%
%       MU(x) = 1 - (2*Re(Lc(x) - Lr(x))) / (Lc(x)^2 - Lr(x)^2) [case C]
%       MU(x) = 1 - 2 / (Lmax(x) + Lmin(x))                     [case R]
%
%       where Lc and Lr refer to the complex and real eigenvalues,
%       respectively, at some point.  Case R means all eigenvalues at x are
%       real, while case C means that two eigenvalues at x form a complex
%       conjugate pair.
%
%   * 'alternating'
%
%       Alternating control values [1].  The values are equal midrange
%       value corresponding to a set of gamma percentiles
%       (PrcMidAlternating option) over the spatial domain.
%
%   * 'global-midrange'
%
%       Global midrange value [1].  Constant such that error contracts over
%       all points in the spatial domain.  If the controllability condition
%       [1] is violated at some points, gamma values at those are ignored
%       while determining the control parameter value.
%
%   * [scalar | M-vector | NX x NY x NZ]
%
%       User-specified (constant | alternating | spatially variant) feeback
%       control values.
%
% DEPENDENCIES
%
%   dvf.sizeVf
%   util.parseOptArgs
%
% REFERENCES
%
%   [1] A. Dubey*, A.-S. Iliopoulos*, X. Sun, F.-F. Yin, and L. Ren,
%   "Iterative inversion of deformation vector fields with feedback
%   control," Medical Physics, vol. 45, no. 7, pp. 3147-3160, May 2018.
%   DOI: 10.1002/mp.12962.
%
%   [2] A. Dubey, "Symmetric completion of deformable registration via
%   bi-residual inversion," PhD thesis, Duke University, Durham, NC, USA.
%
%
% See also      dvf.jacobian, dvf.eigJacobian, dvf.inversion
%
    
    
    %% PARAMETERS
    
    % spatial domain size and dimensionality
    [szDom, dim] = dvf.sizeVf( Lambda );
        
    % options and default values
    opt.Mask              = true( szDom );
    opt.windowMidrange    = 1;
    opt.prcMidAlternating = [98 50];
    opt.fallbackValue     = 0.8;
    
    % parse optional arguments
    opt = util.parseOptArgs( opt, varargin{:} );
    
    % optional output flags
    flag.P = (nargout > 1);
    
    
    %% INITIALIZATION
    
    % if scalar (isotropic) window size was input, extend it
    if isscalar( opt.windowMidrange )
        opt.windowMidrange = repmat( opt.windowMidrange, [1 dim] );
    end
    
    % expand neighborhood window
    if isvector( opt.windowMidrange )
        opt.windowMidrange = true( opt.windowMidrange );
    end
    
    
    %% FEEDBACK CONTROL PARAMETER VALUES
    
    % calculate feedback control parameter values if NTDC-adaptive control is
    % specified (otherwise, just use the input values as-is)
    
    if ischar( control )                % ======== ADAPTIVE CONTROL
        
        % calculate gamma(x) at all x
        Gamma   = min( real(Lambda) ./ abs(Lambda).^2, [], dim+1 );
        MaskBad = isnan( Gamma ) | isinf( Gamma );
        
        % controllable points mask
        MaskCtrl = (Gamma > 0) & ~MaskBad;
        
        % switch among control schemes
        switch lower(control)
            
          case 'midrange'               % ---- local midrange
            
            Mu = localMidrangeMu( Gamma, MaskCtrl, opt.fallbackValue, ...
                                  opt.windowMidrange );
            
          case 'optimal'                % ---- pointwise optimal
            
            Mu = optimalMu( Lambda, MaskCtrl, opt.fallbackValue, szDom, dim );
            
          case 'alternating'            % ---- global alternating
            
            Mu = 1 - prctile( Gamma(MaskCtrl & opt.Mask), ...
                              100 - opt.prcMidAlternating );
            Mu = shiftdim( Mu(:), -(dim+1) );
            
          case 'global-midrange'        % ---- global midrange
            
            Mu = 1 - min( Gamma(MaskCtrl & opt.Mask) );
            
          otherwise                     % ---- error
            
            error( [mfilename ':InvalidControl'], ...
                   'Invalid control scheme: %s', control );
            
        end  % switch (control scheme)
        
        % limit negative control values to avoid translation component blow-up
        Mu( Mu <= -1 ) = -0.99;
        
    else                                % ======== PREFIXED CONTROL
        
        Mu = control;
        
    end  % if (adaptive/prefixed feedback control parameter)
    
    % control value sanitization (just in case of unexpected numerical issues)
    MaskBad = isnan( Mu ) | isinf( Mu );
    if any( MaskBad, 'all' )
        warning( [mfilename ':BadMuValues'], ...
                 ['NaN or Inf control parameter values found. ' ...
                  'Setting them to %.2f'], opt.fallbackValue );
        Mu(MaskBad) = opt.fallbackValue;
    end
    
    
    %% ERROR PROPAGATION FACTORS
    
    if flag.P
        
        Lambda = reshape( Lambda, [szDom, dim] );
        P      = max( abs( 1 - (1 - Mu) .* Lambda ), [], dim+1 );
        P      = reshape( P, [szDom, size(P,dim+2)] );
        
    end  % if (error propagation factor output?)
    
    
    %% TERMINATION
    
    Mu = squeeze( Mu );
    
    
end



%% ==================================================
%  LOCAL MIDRANGE CONTROL PARAMETER VALUES

function Mu = localMidrangeMu (Gamma, MaskCtrl, val, nhood)
% IN    Gamma           reciprocal spectral gap [NX x NY (x NZ)]
%       MaskCtrl        controllable domain mask[NX x NY (x NZ)] (logical)
%       val             fallback control value  [scalar]
%       nhood           local neighborhood mask [WX x WY (x WZ)] (logical)
% OUT   Mu              local midrange values   [NX x NY (x NZ)]
    
    % set invalid Gamma values to the valid upper bound (1), to avoid
    % messing up local-minimum filter calculations
    Gamma( isnan(Gamma) | isinf(Gamma) ) = 1;
    
    % local-neighborhood Gamma values (IMERODE only supports GPU
    % computations for integer/binary i
    flagGpu = isa( Gamma, 'gpuArray' );
    Gamma   = imerode( gather(Gamma), nhood, 'same' );
    if flagGpu
        Gamma = gpuArray( Gamma );
    end
    
    % local midrange control values
    Mu = 1 - Gamma;
    
    % fallback control value over incontrollable sub-domain
    Mu( ~MaskCtrl ) = val;
    
end



%% ==================================================
%  POINTWISE OPTIMAL CONTROL PARAMETER VALUES

function Mu = optimalMu (Lambda, MaskCtrl, val, szDom, dim)
% IN    Lambda          DVF Jacobian eigvals    [NX x NY (x NZ) x (2|3)]
%       MaskCtrl        controllable domain mask[NX x NY (x NZ)] (logical)
%       val             fallback control value  [scalar]
%       szDom           domain dimensions       [NX NY (NZ)]
%       dim             dimensionality          [2|3]
% OUT   Mu              pointwise optimal values[NX x NY (x NZ)]
    
    % parameters
    tauZero  = 1e-3;
    
    % initialization
    Lambda   = reshape( Lambda, [prod(szDom), dim] );
    MaskCtrl = reshape( MaskCtrl, [prod(szDom), 1] );
    Mu       = zeros( [prod(szDom), 1], 'like', real(Lambda) );
    
    % 2D/3D DVF
    switch dim
        
      case 2                            % ======== 2D (real/conjugate pair)
        
        Mu(MaskCtrl) = 1 - 2 ./ sum( real(Lambda(MaskCtrl,:)), 2 );
        
      case 3                            % ======== 3D (real & complex lambda)
        
        % ----- REAL EIGENVALUES (CASE R)
        
        MaskReal     = (max( abs(imag(Lambda)), [], 2 ) < tauZero) & MaskCtrl;
        Mu(MaskReal) = 1 - 2 ./ (max( real(Lambda(MaskReal,:)), [], 2 ) + ...
                                 min( real(Lambda(MaskReal,:)), [], 2 ));
        
        % ----- COMPLEX EIGENVALUES (CASE C)
        
        % complex-eigenvalues domain mask
        MaskComplex = ~MaskReal & MaskCtrl;
        
        % row-column indices of real and (upper halfspace) complex eigenvalues
        [~, jR] = min( abs(imag(Lambda(MaskComplex,:))), [], 2 );% real col-idx
        [~, jC] = max(     imag(Lambda(MaskComplex,:)) , [], 2 );% complex col-idx
        iCRow   = find( MaskComplex ); % row-idx with complex eigvals
        
        % pointwise complex and real eigenvalue pairs
        LR     = Lambda( sub2ind( size(Lambda), iCRow, jR ) );
        LC     = Lambda( sub2ind( size(Lambda), iCRow, jC ) );
        GammaC = real(LC) ./ abs(LC).^2;
        
        % optimal control values: 2 cases
        MaskC1 = abs( 1 - GammaC .* LR ) <= abs( 1 - GammaC .* LC );
        MuC1   = 1 - GammaC;
        MuC2   = 1 - (2 * real( LC - LR ) ./ (abs(LC) + LR) ./ (abs(LC) - LR));
        
        % sanitize case-C2 values (NaN|Inf can propagate despite masking)
        MuC2( isnan(MuC2) | isinf(MuC2) ) = 0;
        
        % set complex-domain Mu
        Mu(MaskComplex) = (MaskC1 .* MuC1) + (~MaskC1 .* MuC2);
        
      otherwise                         % ======== error
        
        error( [mfilename ':InvalidDimensionality'], ...
               'Only 2D and 3D DVFs are supported' );
        
    end  % switch (dimensionality)
    
    % incontrollable subdomain
    Mu(~MaskCtrl) = val;
    
    % reshape domain size
    Mu = reshape( Mu, szDom );
    
end



%%------------------------------------------------------------
%
% AUTHORS
%
%   Alexandros-Stavros Iliopoulos       ailiop@cs.duke.edu
%   Abhishek Dubey                      abhisdub@cs.duke.edu
%   Xiaobai Sun                         xiaobai@cs.duke.edu
%
% VERSION
%
%   1.0.0 - October 31, 2018
%
% ------------------------------------------------------------
