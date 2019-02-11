function I = imageGridSmooth (sizeI, width, bufferBnd)
%
% IMAGEGRIDSMOOTH - Smooth grid-like image generation
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
%   I = IMAGEGRIDSMOOTH( SIZE )
%   I = IMAGEGRIDSMOOTH( SIZE, WIDTH )
%   I = IMAGEGRIDSMOOTH( SIZE, WIDTH, BUFFER )
%
% INPUT
%
%   SIZE        Size of output image            [D-vector]
%               [D1, D2, ...]
%   WIDTH       Grid spacing (sine half-period) [D-vector|scalar]
%               {default: SIZE/10}
%   BUFFER      Dimension-wise buffer from      [D-vector|scalar]
%               image to boundary
%               {default: 0}
%
% OUTPUT
%
%   I           Smooth-grid image               [D1-by-D2-by-...]
%
% DESCRIPTION
%
%   I = IMAGEGRIDSMOOTH(SIZE,WIDTH,BUFFER) generates an image with a
%   D-dimensional smooth grid-like pattern.  The image is a Kronecker
%   product of 1-dimensional sinusoids with specified half-periods.  The
%   image intensities are normalized to [0,1] range.
%
%   If BUFFER is non-zero for the d-th dimension, then BUFFER(d) points are
%   masked off on both ends of the d-th axis.
%
%
% See also      dvf.imdeform
%
    
    
    %% DEFAULTS
    
    if ~exist( 'width', 'var' ) || isempty( width )
        width = sizeI / 10;
    end
    if~exist( 'bufferBnd', 'var' ) || isempty( bufferBnd )
        bufferBnd = 0;
    end
    
    
    %% INITIALIZATION
    
    % number of dimensions
    nDim = length( sizeI );
    
    % if single width/buffer width was input, use the same for all dimensions
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
    
    
    %% DIMENSION-WISE (1D) IMAGE PROFILE AND BUFFER-MASKS
    
    img1d = arrayfun( @(w,d,sz) ...
                      reshape( sin( pi*w * linspace(0,1,sz) ), ...
                               [ones(1,d-1), sz, 1] ), ...
                      width, (1:nDim), sizeI, ...
                      'UniformOutput', false );
    
    mask1d = arrayfun( @(b,d,sz) ...
                       reshape( [zeros(1,b) ones(1,sz-2*b) zeros(1,b)], ...
                                [ones(1,d-1), sz, 1] ), ...
                       bufferBnd, (1:nDim), sizeI, ...
                       'UniformOutput', false );
    
    
    %% MULTI-DIMENSIONAL IMAGE (WITH BUFFER-ZONES)
    
    % initialize image and buffer-zone mask
    I = img1d{1};
    M = mask1d{1};
    
    % Kronecker product of all remaining 1D profiles
    for d = 2 : nDim
        I = I .* img1d{d};
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
% RELEASE
%
%   1.0.2 - February 11, 2019
%
% CHANGELOG
%
%   1.0.2 (Feb 11, 2019) - Alexandros
%       . changed 1D image profile from binary-valued (sharp edges) to
%         sinusoidal (smooth)
%       . multiplicative composition of 1D image profiles into ND image
%         (used to be additive)
%       . renamed function: imageGrid --> imageGridSmooth
%
%   1.0.0 (Oct 31, 2018) - Alexandros
%       . initial release
%
% ------------------------------------------------------------

