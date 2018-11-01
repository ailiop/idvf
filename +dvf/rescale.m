function Fs = rescale (F, s, interpMethod)
%
% RESCALE - Rescale deformation vector field
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
%   FS = RESCALE( F, S )
%   FS = RESCALE( F, S, INTERP )
%
% INPUT
%
%   F           Deformation vector field        [NX x NY (x NZ) x (2|3)]
%   S           Scale factor or rescaled domain [scalar | SX SY (SZ)]
%               size
%   INTERP      Interpolation method            [string]
%               (see imresize3)
%               {'triangle'}
%
% OUTPUT
%
%   FS          Rescaled DVF                    [RX x RY (x RZ) x (2|3)]
%               (fractional array size -> floor)
%
% DESCRIPTION
%
%   FS = RESCALE(F,S,...) scales the spatial domain over which F is defined
%   by a factor of S.  The spatial dimensions of the array are resized and
%   the (interpolated) vector lengths over the resized domain are scaled
%   accordingly.
%
%   If S is a scalar, it is treated as an isotropic scaling factor, such
%   that the spatial domain size is ceil([NX*S NY*S NZ*S]).  The
%   displacement vectors are scaled by ceil(NX*S)/NX (similarly for the Y
%   and Z dimensions).
%
%   If S is a (2|3)-vector, it is treated as the size of the rescaled
%   spatial domain.  The displacement vectors are scaled by SX/NX
%   (similarly for the Y and Z dimensions).
%
% DEPENDENCIES
%
%   dvf.sizeVf
%
%
% See also      imresize3, dvf.inversion
%
    
    
    %% PARAMETERS
    
    % default values for optional arguments
    if ~exist( 'interpMethod', 'var' ) || isempty( interpMethod )
        interpMethod = 'linear';
    end
    
    
    %% INITIALIZATION
    
    % DVF domain size and dimensionality
    [szDom, dim] = dvf.sizeVf( F );
    if length(szDom) ~= dim
        error( [mfilename ':UnmatchedDimensions'], ...
               'Spatial domain size must match vector dimensionality' );
    end
    
    % if no scaling is needed, just exit
    if isequal( s, 1 ) || isequal( s, szDom )
        Fs = F;
        return
    end
    
    % input: scaling factor or scaled-domain size?
    if isscalar( s )                    % ---- scaling factor
        szDomScl = ceil( szDom * s );
    else                                % ---- scaled-domain size
        szDomScl = s;
    end
    s = szDomScl ./ szDom;
    
    % preallocate space for rescaled DVF array
    Fs = zeros( [szDomScl, dim], 'like', F );
    
    % IMRESIZE (2D/3D) annoyingly only supports cubic interpolation on
    % the GPU
    if isa( F, 'gpuArray' )
        interpMethod = 'cubic';
    end
    
    % IMRESIZE (2D) also annoyingly uses different string IDs for some
    % interpolation methods; make that translation as necessary
    if dim == 2
        interpMethod = renameInterpImresize2( interpMethod );
    end
    
    
    %% RESCALE DVF
    
    switch dim
        
      case 2                            % ==== 2D
        
        s = reshape( s, [1 1 dim] );
        for d = 1 : dim
            Fs(:,:,d) = imresize( F(:,:,d), szDomScl, interpMethod );
        end
        
      case 3                            % ==== 3D
        
        s = reshape( s, [1 1 1 dim] );
        for d = 1 : dim
            Fs(:,:,:,d) = imresize3( F(:,:,:,d), szDomScl, interpMethod );
        end
        
      otherwise                         % ==== error
        
        error( [mfilename ':InvalidDimensionality'], ...
               'Only 2D and 3D DVFs are supported' );
        
    end  % switch (dimensionality)
    
    Fs = Fs .* s;
    
    
end



%% ==================================================
%  IMRESIZE 2D INTERPOLATION METHOD STRING ID TRANSLATION

function minterp2 = renameInterpImresize2 (minterpn)
% IN    minterpn        generic interpolation ID        [string]
% OUT   minterp2        2D imresize interpolation ID    [string]
    
    switch lower(minterpn)
        
      case 'linear'
        minterp2 = 'bilinear';
        
      case 'cubic'
        minterp2 = 'bicubic';
        
      otherwise
        minterp2 = minterpn;
        
    end  % switch
    
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

