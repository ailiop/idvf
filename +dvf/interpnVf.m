function Fq = interpnVf (F, gridQ, mInterp, vExtrap)
%
% INTERPNVF - N-dimensional vector field interpolation
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
%   FQ = INTERPNVF( F, GRIDQ )
%   FQ = INTERPNVF( F, GRIDQ, INTERP )
%   FQ = INTERPNVF( F, GRIDQ, INTERP, EXTRAP )
%
% INPUT
%
%   F           Displacement vector field       [N1 x ... x ND x D]
%   GRIDQ       Query points                    [D-cell(M1x...xMD)]
%   INTERP      Interpolation method            [string]
%               {'linear'}
%   EXTRAP      Extrapolation value             [scalar]
%               {NaN}
%
% OUTPUT
%
%   FQ          Interpolated DVF at query points[M1 x ... x MD x D]
%
% DESCRIPTION
%
%   INTERPNVF is a wrapper over the MATLAB built-in function interpn.  It
%   simply loops over all vector (codomain) dimensions and interpolates the
%   corresponding scalar/component fields via a call to interpn.
%
% DEPENDENCIES
%
%   dvf.sizeVf
%
%
% See also      interpn
%
    
    
    %% PARAMETERS
    
    % default values for optional arguments
    if ~exist( 'mInterp', 'var' ) || isempty( mInterp )
        mInterp = 'linear';
    end
    if ~exist( 'vExtrap', 'var' ) || isempty( vExtrap )
        vExtrap = NaN;
    end
    
    % GPU-supported interpolation methods
    mInterpGpu = {'linear', 'nearest'};
    
    
    %% INITIALIZATION
    
    % make sure the interpolation method is supported if running on the GPU
    if isa( F, 'gpuArray' ) && ~any( strcmpi( mInterp, mInterpGpu ) )
        mInterp = mInterpGpu{1};
    end
    
    % reference domain size and DVF dimensionality
    [szDomR, dim] = dvf.sizeVf( F );
    
    % query domain size
    szDomQ = cellfun( @(x,d) size(x,d), gridQ, num2cell(1:dim) );
    
    % linearize spatial dimensions (for support of any dimensionality)
    F     = reshape( F, [prod(szDomR), dim] );
    gridQ = cellfun( @(x) reshape( x, [], 1 ), gridQ, 'UniformOutput', false );
    
    
    %% DVF INTERPOLATION
    
    % interpolate each dimension/scalar component of the DVF
    Fq = zeros( [prod(szDomQ), dim], 'like', F );
    for d = 1 : dim
        Fq(:,d) = interpn( reshape( F(:,d), szDomR ), gridQ{:}, ...
                           mInterp, vExtrap );
    end  % for (d)
    
    
    %% TERMINATION
    
    % re-fold query spatial domain
    Fq = reshape( Fq, [szDomQ, dim] );
    
    
end



%%------------------------------------------------------------
%
% AUTHORS
%
%   Alexandros-Stavros Iliopoulos       ailiop@cs.duke.edu
%
% RELEASE
%
%   1.0.0 - October 31, 2018
%
% ------------------------------------------------------------


