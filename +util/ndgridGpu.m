function varargout = ndgridGpu (gridVec, flagGpu)
%
% NDGRIDGPU - NDGRID wrapper with size-input and gpuArray output options
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
%   [XX, YY, ...] = NDGRIDCELL( {VX, VY, ...} )
%   [XX, YY, ...] = NDGRIDCELL( SIZE )
%   [XX, YY, ...] = NDGRIDCELL( ..., FLAGGPU )
%
% INPUT
%
%   SIZE        Grid domain size                [NX NY ...]
%   VX|VY|...   Dimension-wise grid vectors     [(NX|NY|...)-vector]
%   FLAGGPU     Generate gpuArray grid?         [boolean]
%               {true, iff input is gpuArray}
%
% OUTPUT
%
%   XX|YY|...   X|Y|...-axis grid coordinate    [NX x NY x ...]
%               arrays
%
% DESCRIPTION
%
%   [XX,YY,...] = NDGRIDCELL({VX,VY,...}) is a wrapper around [XX,YY,...] =
%   ndgrid(VX,VY,...).
%
%   [XX,YY,...] = NDGRIDCELL(SIZE) is a similar wrapper, where the default
%   grid vectors are used: [XX,YY,...] = ndgrid(1:SIZE(1),1:SIZE(2),...).
%
%   [XX,YY,...] = NDGRIDCELL(...,true) returns the output grid coordinate
%   vectors as gpuArray objects.
%
% DEPENDENCIES
%
%   <none>
%
%
% See also      ndgrid
%
    
    
    %% PARAMETERS
    
    % main (CPU) memory data by default
    if ~exist( 'flagGpu', 'var' ) || isempty( flagGpu )
        flagGpu = isa( gridVec, 'gpuArray' ) ...
                  || (iscell( gridVec ) && isa( gridVec{1}, 'gpuArray' ));
    end
    
    
    %% INITIALIZATION
    
    % initialize varargout cell
    varargout = cell( size(gridVec) );
    
    % conditional gpuArray function helpers for convenience
    if flagGpu
        hGpu = @gpuArray;
    else
        hGpu = @(x) x;
    end
    
    
    %% ND-GRID
    
    % prepare grid vectors
    if iscell( gridVec )                % ---- grid vectors input
        gridVec = cellfun( hGpu, gridVec, 'UniformOutput', false );
    else                                % ---- grid size input
        gridVec = arrayfun( @(n) hGpu(1):hGpu(n), gridVec, ...
                            'UniformOutput', false );
    end  % if (grid vector/size input?)
    
    % generate ND-grid
    [varargout{:}] = ndgrid( gridVec{:} );
    
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
