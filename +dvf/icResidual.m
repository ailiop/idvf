function R = icResidual (F, G, interpMethod, extrapVal, xxF, yyF, zzF)
%
% ICRESIDUAL - Inverse consistency residual field calculation
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
%   R = ICRESIDUAL( F, G )
%   R = ICRESIDUAL( F, G, INTERP )
%   R = ICRESIDUAL( F, G, INTERP, EXTRAP )
%   R = ICRESIDUAL( F, G, INTERP, EXTRAP, XXF, YYF )
%   R = ICRESIDUAL( F, G, INTERP, EXTRAP, XXF, YYF, ZZF )
%
% INPUT
%
%   F           Forward DVF                     [NX x NY (x NZ) x (2|3)]
%   G           Inverse DVF                     [NX x NY (x NZ) x (2|3)]
%   INTERP      Interpolation method            [string]
%               (see interpn)
%               {'linear'}
%   EXTRAP      Extrapolation value (for out-   [scalar]
%               of-boundary interpolation)
%               (see interpn)
%               {NaN}
%   XXF|YYF|ZZF Forward-DVF displaced           [NX x NY (x NZ)]
%               coordinates of IC residual
%               spatial domain
%               {dvf.coorDisplace(F)}
%
% OUTPUT
%
%   R           Inverse consistency residual    [NX x NY (x NZ) x (2|3)]
%
% DESCRIPTION
%
%   R = ICRESIDUAL(F,G) computes the inverse consistency (IC) residual
%   field [1-3]: 
%
%       R(x) = F(x) + G(x + F(x)) .
%
%   R = ICRESIDUAL(F,G,INTERP) and 
%   R = ICRESIDUAL(F,G,INTERP,EXTRAP) allow the user to specify the
%   interpolation method and default extrapolation value to be used when
%   computing the G(x + F(x)) values.
%
%   R = ICRESIDUAL(F,G,INTERP,EXTRAP,XXF,YYF,ZZF) evaluates G at the
%   user-supplied 3D coordinates in XXF, YYF, and ZZF, as per
%
%       R = F + interpn( G, XXF, YYF, ZZF, ... )
%
%   for each vector-field (codomain) dimension, and similarly for 2D DVFs.
%
% DEPENDENCIES
%
%   dvf.coorDisplace
%   dvf.interpnVf
%
% REFERENCES
%
%   [1] A. Dubey*, A.-S. Iliopoulos*, X. Sun, F.-F. Yin, and L. Ren,
%   "Iterative inversion of deformation vector fields with feedback
%   control," Medical Physics, vol. 45, no. 7, pp. 3147-3160, 2018.
%   [DOI:10.1002/mp.12962]
%
%   [2] A. Dubey, "Symmetric completion of deformable registration via
%   bi-residual inversion," PhD thesis, Duke University, Durham, NC, USA.
%
%   [3] G. E. Christensen and H. J. Johnson, "Consistent image
%   registration," IEEE Transactions on Medical Imaging, vol. 20, no. 7,
%   pp. 568-582, 2001.
%
%
% See also      dvf.inversion, dvf.coorDisplace, dvf.interpnVf
%
    
   
    %% PARAMETERS
    
    % DVF dimensionality
    [~, dim] = dvf.sizeVf( F );
    
    % default values for optional arguments
    if ~exist( 'interpMethod', 'var' ) || isempty( interpMethod )
        interpMethod = 'linear';
    end
    if ~exist( 'extrapVal', 'var' ) || isempty( extrapVal )
        extrapVal = NaN;
    end
    if (dim == 2) ...
             && (~exist( 'xxF', 'var' ) || isempty( xxF ) || ...
                 ~exist( 'yyF', 'var' ) || isempty( yyF ))
        [xxF, yyF] = dvf.coorDisplace( F );
    end
    if (dim == 3) ...
             && (~exist( 'xxF', 'var' ) || isempty( xxF ) || ...
                 ~exist( 'yyF', 'var' ) || isempty( yyF ) || ...
                 ~exist( 'zzF', 'var' ) || isempty( zzF ))
        [xxF, yyF, zzF] = dvf.coorDisplace( F );
    end
    
    
    %% IC RESIDUAL COMPUTATION
    
    switch dim
        
      case 2                            % -------- 2D
        
        R = F + dvf.interpnVf( G, {xxF, yyF}, interpMethod, extrapVal );
        
      case 3                            % -------- 3D
        
        R = F + dvf.interpnVf( G, {xxF, yyF, zzF}, interpMethod, extrapVal );
        
      otherwise                         % -------- error
        
        error( [mfilename ':InvalidDimensionality'], ...
               'Only 2D and 3D DVFs are supported' );
        
    end  % switch (dimensionality)
    
    
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
