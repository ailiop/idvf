function Mask = maskDomain (F)
%
% MASKDOMAIN - Valid displacement domain mask
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
%   MASK = MASKDOMAIN( F )
%
% INPUT
%
%   F           Displacement vector field       [NX x NY (x NZ) x (2|3)]
%
% OUTPUT
%
%   MASK        Valid displacement domain mask  [NX x NY (x NZ)]
%
% DESCRIPTION
%
%   MASK = MASKDOMAIN(F) returns a mask over the "valid" displacement
%   domain, defined as the common or non-occluded intersection of the input
%   and deformed spatial domain [1]:
%
%       Mask = Domain - (MaskOut \union MaskNoOvlp)
%
%   where       MaskOut    = {x | x \in Domain  &  DVF(x) \notin Domain} ,
%               MaskNoOvlp = Domain - {DVF(x) | x \in Domain} ,
%
%   and Domain is the spatial domain of the input DVF.
%
% DEPENDENCIES
%
%   dvf.coorDisplace
%   util.ndgridGpu
%
% REFERENCES
%
%   [1] A. Dubey*, A.-S. Iliopoulos*, X. Sun, F.-F. Yin, and L. Ren,
%   "Iterative inversion of deformation vector fields with feedback
%   control," Medical Physics, vol. 45, no. 7, pp. 3147-3160, May 2018.
%   DOI: 10.1002/mp.12962.
%
%
% See also      dvf.inversion, dvf.coorDisplace, imclose
%
    
    
    %% INITIALIZATION
    
    % DVF domain size and dimensionality
    [szDom, dim] = dvf.sizeVf( F );
    
    % dimensionality check
    if (dim ~= 2) && (dim ~= 3)
        error( [mfilename ':InvalidDimensionality'], ...
               'Only 2D and 3D vector fields are supported' );
    end
    
    % max displacement in each dimension [1 x D]
    maxDisp = ceil( reshape( max( max( abs(F), [], 1:dim ), 0 ), [1 dim] ) );
    
    
    %% DOMAIN MASK
    
    % deformed-grid coordinates
    gridDfm      = cell( 1, dim );
    [gridDfm{:}] = dvf.coorDisplace( F );
    
    % within-bounds domain
    MaskIn = true( szDom );
    for d = 1 : dim
        MaskIn = MaskIn & (gridDfm{d} >= 1) & (gridDfm{d} <= szDom(d));
    end
    
    % deformed domain (equivalent to nearest-neighbor scattered interpolation
    % from deformed coordinates onto padded reference grid)
    MaskDfm = false( szDom + 2*maxDisp + 1 );
    gridDfm = cellfun( @(x,d) round(x) + d, gridDfm, num2cell(maxDisp), ...
                       'UniformOutput', false );
    idxDfm  = sub2ind( size(MaskDfm), gridDfm{:} );
    MaskDfm(idxDfm) = true;
    
    % morphological image closure to fill "holes" and "canals"
    MaskDfm = imclose( MaskDfm, true( 2*maxDisp+1 ) );
    
    % depad mask (equivalent to intersection with the reference domain)
    depad   = arrayfun( @(sz,d) (1 : sz) + d, szDom, maxDisp, ...
                        'UniformOutput', false );
    MaskDfm = MaskDfm(depad{:});
    
    % within-bounds and deformed domain intersection
    Mask = MaskIn & MaskDfm;
    
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
