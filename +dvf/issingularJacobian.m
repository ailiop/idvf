function MaskSingular = issingularJacobian (J, varargin)
%
% ISSINGULARJACOBIAN - Mask of singular Jacobian fields
%   
% SYNTAX
%
%   MASKSNG = ISSINGULARJACOBIAN( J )
%   MASKSNG = ISSINGULARJACOBIAN( J, 'Name', Value, ... )
%
% INPUT
%
%   J           Deformation Jacobians           [NX x NY (x NZ) x (2|3)]
%
% OUTPUT
%
%   MASKSNG     Numerically singular Jacobians  [NX x NY (x NZ)] (logical)
%               mask
%
% OPTIONS
%
%   'Threshold'         [scalar]                        {1e-6}
%
%       Numerical zero threshold for singularity test via condition
%       number or determinant of the Jacobians.  Positive scalar.
%
%   'ExplicitDet'       [boolean]                       {false}
%
%       false => check the condition number of each Jacobian, implemented
%       via looping COND over the J array slices.
%
%       true => check the determinant of each Jacobian, implemented via
%       explicit 2D/3D determinant calculation over the J array slices.
%       WARNING: this option is faster in the case of a large number of
%       Jacobians, but it is numerically unstable.  (Only works with 2D/3D
%       Jacobians.)
%
%   'Mask'              [NX x NY (x NZ)]                {true(NX,NY(,NZ))}
%
%       Logical array indicating the spatial sub-domain of interest.
%       Singularity tests are only done over points where the input Mask is
%       true; the output MASKSNG is automatically false over the rest of
%       the domain.
%
% DESCRIPTION
%
%   MASKSNG = ISSINGULARJACOBIAN(J,...) returns a mask which is set over
%   all points in the spatial domain where the deformation Jacobian is
%   numerically singular.  The DVF is non-invertible over the sub-domain
%   indicated by MASKSNG.
%
% DEPENDENCIES
%
%   util.parseOptArgs
%
%
% See also      dvf.inversion
%
    
    
    %% PARAMETERS
    
    % domain size and dimensionality
    szDom = size( J );
    dim   = szDom(end);
    szDom = szDom(1:end-2);
    
    % options and default values
    opt.threshold   = 1e-6;
    opt.explicitDet = false;
    opt.Mask        = true( szDom );
    
    % parse optional arguments
    opt = util.parseOptArgs( opt, varargin{:} );
    
    
    %% INITIALIZATION
    
    % linearize and mask spatial domain
    JMask = reshape( J, [prod(szDom), dim, dim] );
    JMask = JMask(opt.Mask,:,:);  % [NxDxD]
    
    
    %% SINGULARITY TESTS
    
    if opt.explicitDet                  % ======== determinant
        
        detJMask = detJacobian( JMask, dim );
        Mask     = (abs(detJMask) <= opt.threshold);
        
    else                                % ======== condition number
        
        JMask     = permute( JMask, [2 3 1] );  % [DxDxN]
        condJMask = zeros( [size(JMask,3), 1], 'like', J );
        for i = 1 : length(condJMask)
            condJMask(i) = cond( JMask(:,:,i) );
        end
        Mask = (condJMask >= 1/opt.threshold);
        
    end  % if (determinant / condition number calculation)
    
    
    %% TERMINATION
    
    % inject subdomain mask to whole-domain mask
    MaskSingular           = false( [prod(szDom), 1] );
    MaskSingular(opt.Mask) = Mask;
    MaskSingular           = reshape( MaskSingular, [szDom 1] );
    
    
end



%% ==================================================
%  2D/3D DETERMINANT (EXPLICIT CALCULATION)

function detJ = detJacobian (J)
% IN    J               Jacobian matrices       [N x D x D]
% OUT   detJ            Jacobian determinants   [N x 1]
    
    switch size(J,3)
        
      case 2                            % ==== 2D
        
        detJ = J(:,1,1) * J(:,2,2) - J(:,1,2) * J(:,2,1);
        
      case 3                            % ==== 3D
        
        detJ = J(:,1,1) * (J(:,2,2) * J(:,3,3) - ...
                           J(:,2,3) * J(:,3,2)) - ...
               J(:,1,2) * (J(:,2,1) * J(:,3,3) - ...
                           J(:,2,3) * J(:,3,1)) + ...
               J(:,1,3) * (J(:,2,1) * J(:,3,2) - ...
                           J(:,2,2) * J(:,3,1));
        
      otherwise                         % ==== error
        
        error( [mfilename ':InvalidDimensionality'], ...
               ['Explicit determinant calculations only supported with' ...
                ' 2D/3D Jacobians'] );
        
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
%   1.0.1 - November 7, 2018
%
% ------------------------------------------------------------
