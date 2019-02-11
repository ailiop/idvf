function IDfm = imdeform (IRef, DVF, mapping, interpMethod, extrapVal)
%
% IMDEFORM - 2D/3D image deformation
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
%   IDFM = IMDEFORM( IREF, DVF )
%   IDFM = IMDEFORM( IREF, DVF, MAPPING )
%   IDFM = IMDEFORM( IREF, DVF, MAPPING, MINTERP )
%   IDFM = IMDEFORM( IREF, DVF, MAPPING, MINTERP, VEXTRAP )
%
% INPUT
%
%   IREF        Reference (to-be-deformed) image        [NX x NY x NZ]
%   DVF         Deformation vector field                [NX x NY x NZ x 3]
%               (displacement measured in pixels)
%   MAPPING     Deformation mapping                     [string]
%               * 'bwd'|'backward'
%               * 'fwd'|'forward'
%               {default: 'bwd'}
%   MINTERP     Interpolation method                    [string]
%               (see interpn and griddata)
%               {default: 'linear'}
%   VEXTRAP     Extrapolation value                     [scalar]
%               {default: 0}
%
% OUTPUT
%
%   IDFM        Deformed image                          [NX x NY x NZ]
%
% NOTE
%
%   The image/vector field array sizes in the documentation refer to the
%   3D case.  In the 2D case, the sizes are amended accordingly:
%       - [NX x NY x NZ x 3] --> [NX x NY x 2]
%       - [NX x NY x NZ]     --> [NX x NY]
%
% DESCRIPTION
%
%   IDFM = IMDEFORM(IREF,DVF) or 
%   IDFM = IMDEFORM(IREF,DVF,'bwd') deforms the input reference 3D image,
%   such that
%
%               IDFM(x) = IREF(x + DVF(x)) ,
%
%   for all x in the deformed-image domain. The above model is referred to
%   as the backward-mapping deformation model.
%
%   IDFM = IMDEFORM(IREF,DVF,'fwd') applies the forward-mapping deformation
%   model, instead:
%
%               IDFM(x + DVF(x)) = IREF(x) ,
%
%   for all x in the reference-image domain. This is much slower,
%   especially for large images, as it is based on scattered data
%   interpolation.  Instead, the backward-mapping model is based on gridded
%   data interpolation.
%
%   IDFM = IMDEFORM(IREF,DVF,...,MINTERP) also specifies the interpolation
%   method to be used when regridding the deformed image onto the regular
%   Cartesian grid of the reference image domain.
%
%   IDFM = IMDEFORM(IREF,DVF,...,VEXTRAP) also specifies the extrapolaiton
%   value to be used when interpolating outside the bounding box of the
%   reference grid.
%
% DEPENDENCIES
%
%   dvf.sizeVf
%   dvf.coorDisplace
%   util.ndgridGpu
%
%
% See also      interpn, griddata, dvf.coorDisplace
%
    
    
    %% PARAMETERS
    
    % deformation mapping
    if ~exist( 'mapping', 'var' ) || isempty( mapping )
        mapping = 'bwd';
    end
    
    % interpolation method
    if ~exist( 'interpMethod', 'var' ) || isempty( interpMethod )
        interpMethod = 'linear';
    end
    
    % extrapolation method
    if ~exist( 'extrapVal', 'var' ) || isempty( extrapVal )
        extrapVal = 0;
    end
    
    
    %% DIMENSIONALITY
    
    % spatial domain size and dimensionality
    [szDom, dim] = dvf.sizeVf( DVF );
    
    % ensure 2D/3D case (the code actually supports ND image deformation)
    if ~( ((dim == 3) && (length(szDom) == 3)) || ...
          ((dim == 2) && (length(szDom) == 2)) )
        error( [mfilename ':InvalidDimensionality'], ...
               'Only 2D/3D deformation is supported' );
    end
    
    
    %% INITIALIZATION
    
    % reference grid coordinates
    gridRef      = cell( 1, dim );
    [gridRef{:}] = util.ndgridGpu( szDom, isa( DVF, 'gpuArray' ) );
    
    % deformed grid coordinates
    gridDfm      = cell( 1, dim );
    [gridDfm{:}] = dvf.coorDisplace( DVF, gridRef{:} );
    
    
    %% IMAGE DEFORMATION
    
    % deform the reference image
    switch lower( mapping )
        
      case {'bwd','backward'}   % ----- BACKWARD MAPPING
        
        IDfm = interpn( IRef, gridDfm{:}, interpMethod, NaN );
        
      case {'fwd','forward'}    % ----- FORWARD MAPPING
        
        gridDfm = cellfun( @(x) x(:), gridDfm, 'UniformOutput', false );
        IDfm    = griddata( gridDfm{:}, IRef(:), gridRef{:}, interpMethod );
        
      otherwise                 % ----- ERROR (UNKNOWN MAPPING)
        
        error( 'imdeform:UnknownMapping', ...
               'MAPPING must be [''bwd''|''fwd'']' );
        
    end  % switch (deformation mapping)
    
    % set any out-of-domain (NaN) values to 0
    IDfm( isnan(IDfm) ) = extrapVal;
    
    
end




%%------------------------------------------------------------
%
% AUTHORS
%
%   Alexandros-Stavros Iliopoulos       ailiop@cs.duke.edu
%
% RELEASE
%
%   1.0.2 - Ferbuary 11, 2019
%
% CHANGELOG
%
%   1.0.2 (Feb 11, 2019) - Alexandros
%       . documentation clean-up
%
%   1.0.0 (Oct 31, 2018) - Alexandros
%       . initial release
%
% ------------------------------------------------------------
