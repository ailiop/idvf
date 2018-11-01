function [xxd, yyd, zzd] = coorDisplace (Dvf, xx, yy, zz)
%
% COORDISPLACE - Calculate displaced coordinates
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
%   [XXD, YYD] = COORDISPLACE( DVF )
%   [XXD, YYD] = COORDISPLACE( DVF, XX, YY )
%   [XXD, YYD, ZZD] = COORDISPLACE( DVF )
%   [XXD, YYD, ZZD] = COORDISPLACE( DVF, XX, YY, ZZ )
%
% INPUT
%
%   DVF         Displacement vector field       [NX x ... x (2|3)]
%   XX|YY|ZZ    Reference coordinates           [NX x ...]
%
% OUTPUT
%
%   XXD|YYD|ZZD Displaced coordinates           [NX x ...]
%
% DESCRIPTION
%
%   [XXD,YYD,ZZD] = COORDISPLACE(DVF,XX,YY,ZZ) is a simple
%   wrapper around 3 coordinate displacement calls:
%
%       XXD = XX + DVF(:,1) ,
%
%   and similarly along the Y-axis (as well as Z-axis in the 3D case).
%
%   The spatial-domain array dimensions of the input DVF (i.e. all but the
%   trailing dimension) may be in any shape---linearized, 3D ndgrid, or
%   anything else.  The output coordinate arrays will be shaped as the
%   input DVF.  (Note that this array-shape flexibility is not evident in
%   the above coordinates-displacement statement.)
%
%   [XXD,YYD,ZZD] = COORDISPLACE(DVF) applies the displacement
%   field to the default grid:
%
%       [XX,YY,ZZ] = ndgrid( 1:size(DVF,1), 1:size(DVF,2), 1:size(DVF,3) )
%
% DEPENDENCIES
%
%   dvf.sizeVf
%   util.ndgridGpu
%
%
% See also      ndgrid, interpn
%
    
    
    %% DEFAULT GRID?
    
    [szDom, dim] = dvf.sizeVf( Dvf );
    
    if (dim == 2) ...
             && (~exist( 'xx', 'var' ) || isempty( xx ) || ...
                 ~exist( 'yy', 'var' ) || isempty( yy ))
        [xx,yy] = util.ndgridGpu( szDom, isa( Dvf, 'gpuArray' ) );
    end
    if (dim == 3) ...
             && (~exist( 'xx', 'var' ) || isempty( xx ) || ...
                 ~exist( 'yy', 'var' ) || isempty( yy ) || ...
                 ~exist( 'zz', 'var' ) || isempty( zz ))
        [xx,yy,zz] = util.ndgridGpu( szDom, isa( Dvf, 'gpuArray' ) );
    end
    
    
    %% INITIALIZATION
    
    % linearize spatial domain dimensions
    Dvf = reshape( Dvf, [prod(szDom), dim] );
    xx  = reshape( xx , [prod(szDom), 1  ] );
    yy  = reshape( yy , [prod(szDom), 1  ] );
    if dim == 3
        zz = reshape( zz, [prod(szDom), 1  ] );
    end
    
    
    %% DEFORMED GRID
    
    xxd = xx + Dvf(:,1);
    yyd = yy + Dvf(:,2);
    if dim == 3
        zzd = zz + Dvf(:,3);
    end
    
    
    %% TERMINATION
    
    % reshape spatial domain to match input DVF
    xxd = reshape( xxd, [szDom, 1] );
    yyd = reshape( yyd, [szDom, 1] );
    if dim == 3
        zzd = reshape( zzd, [szDom, 1] );
    end
    
    
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

