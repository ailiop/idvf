function [sizeV, dimV] = sizeVf (V)
%
% DIMENSIONS - Vector-field array spatial domain size and dimensionality
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
%   [SIZE, DIM] = DIMENSIONS( VF )
%
% INPUT
%
%   VF          Vector field                    [N1 x ... x NK x DIM]
%
% OUTPUT
%
%   SIZE        Vector field spatial domain size[N1 ... NK]
%   DIM         Vector field dimensionality     [DIM]
%
% DESCRIPTION
%
%   DIMENSIONS is a simple utility function for getting the size/shape of
%   the spatial domain of a vector-field array and its dimensionality.  It
%   simply splits the output of the (K+1)-vector size(VF) to the K-vector
%   head (SIZE) and the last-element tail (DIM).
%
% DEPENDENCIES
%
%   <none>
%
%
% See also      dvf.coorDisplace, dvf.feedbackControlVal
%
    
    
    sizeV = size( V );
    dimV  = sizeV(end);
    sizeV = sizeV(1:end-1);
    
    
end



%%------------------------------------------------------------
%
% AUTHORS
%
%   Alexandros-Stavros Iliopoulos       ailiop@cs.duke.edu
%
% RELEASE
%
%   1.0.0 - October 31, 2018
%
% ------------------------------------------------------------

