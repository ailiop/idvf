function Fb = changeUnit (Fa, atob)
%
% CHANGEUNIT - Change vector-field unit of measurement
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
%   FB = CHANGEUNIT( FA, ATOB )
%
% INPUT
%
%   FA          DVF in A units                  [NX x ... x NK x D]
%   ATOB        A-to-B conversion scale         [scalar | D-vector]
%
% OUTPUT
%
%   FB          DVF in B units                  [NX x ... x NK x D]
%
% DESCRIPTION
%
%   FB = CHANGEUNIT(FA,ATOB) is a simple utility function for changing the
%   input vector-field codomain unit of measurement without rescaling the
%   domain (see dvf.rescale for that).  ATOB is the conversion factor(s);
%   it may a scalar (isotropic conversion) or a D-vector (anisotropic
%   conversion).
%
%   CHANGEUNIT may be useful for inverse consistency/error evaluation: most
%   functions in the +dvf package expect that input DVFs are measured in
%   pixels/voxels, but it may be desirable to measure post-inversion error
%   in some other unit (e.g., mm).
%
% DEPENDENCIES
%
%   dvf.sizeVf
%
%
% See also      dvf.rescale, dvf.icResidual
%
    
    
    % DVF domain size and dimensionality
    [szDom, dim] = dvf.sizeVf( Fa );
    
    % scale DVF codomain units
    Fb = reshape( Fa, [prod(szDom), dim] ) .* reshape( atob, 1, [] );
    Fb = reshape( Fb, [szDom, dim] );
    
    
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
