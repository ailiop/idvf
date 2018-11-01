function [CtrlIdx, Rho, Det, Lambda] = ntdcMeasures (F)
%
% NTDCMEASURES - Spectral characterization measures of deformation
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
%   [CTRLIDX, RHO, DET] = NTDCMEASURES( LAMBDA )
%   [CTRLIDX, RHO, DET, LAMBDA] = NTDCMEASURES( DVF )
%
% INPUT
%
%   LAMBDA      Deformation Jacobian eigenvalues[NX x NY (x NZ) x (2|3)] (complex)
%   DVF         Displacement vector field       [NX x NY (x NZ) x (2|3)] (real)
%
% OUTPUT
%
%   CTRLIDX     Algebraic control index         [NX x NY (x NZ)]
%   RHO         NTDC spectral radius            [NX x NY (x NZ)]
%   DET         Deformation determinant         [NX x NY (x NZ)]
%   LAMBDA      Deformation Jacobian eigenvalues[NX x NY (x NZ) x (2|3)] (complex)
%
% DESCRIPTION
%
%   [CTRLIDX,RHO,DET] = NTDCMEASURES(LAMBDA) processes the eigenvalue field
%   of a deformation Jacobian calculates three spectral measure maps that
%   characterize the non-translational displacement component (NTDC) of the
%   underlying deformation [1].
%
%   [CTRLIDX,RHO,DET,LAMBDA] = NTDCMEASURES(DVF) first computes the
%   deformation Jacobian and its eigevalues using the input DVF, and then
%   calculates the three spectral maps.
%
%   *NOTE* The distinction between the two is made solely by whether or
%   not the input array is of 'real' or 'complex' type.  If the input is
%   to be eigenvalues that happen to be strictly real, it should be
%   passed as complex(LAMBDA).
%   
% MEASURES
%
%   * Algebraic control index
%                                                     Re(lambda_j(x))
%       a(x) = 1 - 2*gamma(x),  where gamma(x) = min ----------------- .
%                                                 j   |lambda_j(x)|^2
%       
%       If a(x) >= 1, then there exists no feasible control parameter value
%       such that the DVF inversion iteration will converge (around x).
%   
%       If a(x) > 0 <==> rho(x) >= 1, the NTDC is non-small and active
%       feedback control is necessary in the inversion iteration.
%
%   * NTDC spectral radius
%
%       This is the spectral radius of the *displacement* Jacobian, whose
%       eigenvalues (lambda') are related to those of the transformation
%       Jacobian (lambda) by: lambda' = lambda - 1.  If the spectral radius
%       is greater than 1, the NTDC is considered non-small (see above).
%
%   * Deformation determinant
%
%       det(J(x)) = prod_j(lambda_j(x)).  
%
%       The deformation Jacobian shows whether the deformed volume is
%       locally expanding (det(J(x)) > 1) or contracting (det(J(x)) < 1).
%
%       The DVF is non-invertible in regions where det(J(x)) = 0.
%
%       Negative values (det(J(x)) < 0) indicate that the volume is locally
%       turned inside out.  Such regions are typically artifacts of a
%       deformable registration algorithm (and often coincide with regions
%       where the controllability condition is violated).
%
% DEPENDENCIES
%
%   dvf.jacobian
%   dvf.sizeVf
%   dvf.eigJacobian
%
% REFERENCES
%
%   [1] A. Dubey*, A.-S. Iliopoulos*, X. Sun, F.-F. Yin, and L. Ren,
%   "Iterative inversion of deformation vector fields with feedback
%   control," Medical Physics, vol. 45, no. 7, pp. 3147-3160, May 2018.
%   DOI: 10.1002/mp.12962.
%
%
% See also      dvf.jacobian, dvf.eigJacobian, dvf.feedbackControlVal
%
    
    
    %% INITIALIZATION
    
    % DVF dimensionality
    [szDom, ~] = dvf.sizeVf( F );
    
    % optional output flags
    flagRho = (nargout > 1);
    flagDet = (nargout > 2);
    
    
    %% SPECTRUM COMPUTATIONS
    
    if isreal( F )                      % ==== DVF input
        Lambda = dvf.eigJacobian( dvf.jacobian( F ) );
    else                                % ==== eigenvalues input
        Lambda = F;
    end  % if (DVF/eigenvalues input?)
    
    
    %% SPECTRAL MEASURES
    
    % algebraic control index
    CtrlIdx = 1 - 2*min( real(Lambda) ./ abs(Lambda).^2, [], length(szDom)+1 );
    
    % NTDC spectral radius
    if flagRho
        Rho = max( abs(Lambda - 1), [], length(szDom)+1 );
    end
    
    % deformation determinant (the determinant is guaranteed to be real; the
    % call to REAL below is to avoid potential issues due to numerical
    % errors)
    if flagDet
        Det = real( prod( Lambda, length(szDom)+1 ) );
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

