function Lambda = eigJacobian (J, Mask)
%
% EIGJACOBIAN - Jacobian field eigenvalues
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
%   LAMBDA = EIGJACOBIAN( J )
%   LAMBDA = EIGJACOBIAN( J, MASK )
%
% INPUT
%
%   J           Jacobian tensor field           [N1 x ... x NK x D x D]
%   MASK        Spatial domain mask             [N1 x ... x NK]
%               {true(N1,...,NK)}
%
% OUTPUT
%
%   LAMBDA      Jacobian eigenvalues            [N1 x ... x NK x D]
%
% DESCRIPTION
%
%   LAMBDA = EIGJACOBIAN(J) computes the eigenvalues of the input Jacobian
%   field at each location of its spatial domain.
%
%   LAMBDA = EIGJACOBIAN(J,MASK) only computes the eigenvalues at those
%   locations specified by MASK.  The output array contains NaN for all
%   other locations.
%
% DEPENDENCIES
%
%   <none>
%
%
% See also      dvf.jacobian, dvf.feedbackControlVal, dvf.ntdcMeasures
%
    
    
    %% PARAMETERS
    
    % domain size and dimensionality
    szDom = size( J );
    dim   = szDom(end);
    szDom = szDom(1:end-2);
    
    % spatial mask covers the whole domain by default
    if ~exist( 'Mask', 'var' ) || isempty( Mask )
        Mask = true( szDom );
    end
    
    
    %% INITIALIZATION
    
    % linearize and mask spatial domain
    JMask = reshape( J, [prod(szDom), dim, dim] );
    JMask = JMask(Mask,:,:);
    JMask = permute( JMask, [2 3 1] ); % [DxDxN]
    
    
    %% EIGENVALUE CALCULATION
    
    LambdaMask = zeros( [nnz(Mask), dim], 'like', J );
    for i = 1 : nnz(Mask)
        LambdaMask(i,:) = eig( JMask(:,:,i) );
    end
    
    
    %% TERMINATION
    
    % put masked-domain eigenvalues to full-domain array
    Lambda         = nan( [prod(szDom), dim], 'like', J );
    Lambda(Mask,:) = LambdaMask;
    Lambda         = reshape( Lambda, [szDom, dim] );
    
    
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
