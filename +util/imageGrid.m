function I = imageGrid (sizeI, width, bufferBnd)
%
% IMAGEGRID - Grid-like 2D image generation
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
%   I = IMAGEGRID( SIZE )
%   I = IMAGEGRID( SIZE, WIDTH )
%   I = IMAGEGRID( SIZE, WIDTH, BUFFER )
%
% INPUT
%
%   SIZE        Size of output image            [D-vector]
%               [D1, D2, ...]
%   WIDTH       Checkered-pattern strips width  [D-vector|scalar]
%               {default: SIZE/21}
%   BUFFER      Dimension-wise buffer from      [D-vector|scalar]
%               checkered-pattern to boundary
%               {default: 0}
%
% OUTPUT
%
%   I           Checkered-pattern image         [D1-by-D2-by-...]
%
% DESCRIPTION
%
%   I = IMAGEGRID(SIZE,WIDTH,BUFFER) generates an image with a
%   D-dimensional checkered pattern as per dimension-wise specified strip
%   widths. The image intensities are normalized such that they are K/D,
%   where K is the number of coinciding strips (for each pixel).
%
%   If BUFFER is non-zero for the d-th dimension, then BUFFER(d) points are
%   masked off on both ends of the d-th axis.
%
%
% See also      dvf.imdeform
%
    
    
    %% DEFAULTS
    
    if ~exist( 'width', 'var' ) || isempty( width )
        width = sizeI / 21;
    end
    if~exist( 'bufferBnd', 'var' ) || isempty( bufferBnd )
        bufferBnd = 0;
    end
    
    
    %% INITIALIZATION
    
    % number of dimensions
    nDim = length( sizeI );
    
    % if single strip/buffer width was input, use the same for all
    % dimensions
    if isscalar( width )
        width = repmat( width, [nDim, 1] );
    end
    if isscalar( bufferBnd )
        bufferBnd = repmat( bufferBnd, [nDim, 1] );
    end
    
    % make sure all inputs are row vectors
    sizeI     = reshape( sizeI    , [1, nDim] );
    width     = reshape( width    , [1, nDim] );
    bufferBnd = reshape( bufferBnd, [1, nDim] );
    
    
    %% DIMENSION-WISE (1D) CHECKERED PATTERNS AND BUFFER-MASKS
    
    chk1d = arrayfun( @(w,d,sz) ...
                      reshape( (mod( 0:(sz-1), 2*w ) < w), ...
                               [ones(1,d-1), sz, 1] ), ...
                      width, (1:nDim), sizeI, ...
                      'UniformOutput', false );
    
    mask1d = arrayfun( @(b,d,sz) ...
                       reshape( [zeros(1,b) ones(1,sz-2*b) zeros(1,b)], ...
                                [ones(1,d-1), sz, 1] ), ...
                       bufferBnd, (1:nDim), sizeI, ...
                       'UniformOutput', false );
    
    
    %% CHECKERED-PATTERN IMAGE (WITH BUFFER-ZONES)
    
    % initialize checkered-pattern image and buffer-zone mask
    I = chk1d{1};
    M = mask1d{1};
    
    % Kronecker summation of all remaining 1D patterns
    for d = 2 : nDim
        I = I  + chk1d{d};
        M = M .* mask1d{d};
    end
    
    % mask-out buffer region points
    I = I .* M;
    
    
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

