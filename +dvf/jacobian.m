function J = jacobian (F)
%
% JACOBIAN - 2D/3D deformation transformation Jacobian tensor field
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
%   J = JACOBIAN( F )
%
% INPUT
%
%   F           Displacement vector field       [N1 x N2 (x N3) x (2|3)]
%
% OUTPUT
%
%   J           Deformation transformation      [NX x NY (x N3) x (2|3) x (2|3)]
%               Jacobian
%
% DESCRIPTION
%
%   J = JACOBIAN(F) computes the Jacobian tensor field of the deformation
%   transformation described by the input DVF.  The transformation Jacobian
%   is related to the displacement Jacobian by [1]
%
%       J = D(F) + I
%
%   where D(F) is the Jacobian (gradient tensor) of the displacement field
%   F, and I is the 2x2 (or 3x3) identiy matrix.
%
%   Calculations are done using central finite differences (see GRADIENT).
%
% DEPENDENCIES
%
%   <none>
%
% REFERENCES
%
%   [1] A. Dubey*, A.-S. Iliopoulos*, X. Sun, F.-F. Yin, and L. Ren,
%   "Iterative inversion of deformation vector fields with feedback
%   control," Medical Physics, vol. 45, no. 7, pp. 3147-3160, May 2018.
%   DOI: 10.1002/mp.12962.
%
% See also      gradient, dvf.inversion, dvf.eigJacobian
%
    
    
    %% INITIALIZATION
    
    % domain size and dimensionality
    [szDom, dim] = dvf.sizeVf( F );
    
    % preallocate space for output Jacobian array
    J = zeros( [szDom dim dim], 'like', F );
    
    
    %% DEFORMATION TRANSFORMATION JACOBIAN
    
    switch dim
        
      case 2                            % -------- 2D
        
        for d = 1 : dim
            [J(:,:,d,2), J(:,:,d,1)] = gradient( F(:,:,d), 1 );
            J(:,:,d,d) = J(:,:,d,d) + 1;
        end
        
      case 3                            % -------- 3D
        
        for d = 1 : dim
            [J(:,:,:,d,2), J(:,:,:,d,1), J(:,:,:,d,3)] = ...
                gradient( F(:,:,:,d), 1 );
            J(:,:,:,d,d) = J(:,:,:,d,d) + 1;
        end
        
      otherwise                         % -------- error
        
        error( [mfilename ':InvalidDimensionality'], ...
               'Only 2D and 3D vector fields are supported.' );
        
    end  % if (DVF dimensionality)
    
    
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

