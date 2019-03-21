%
% SCRIPT: test_idvf.m
%
%
% Test idvf functionality by running the included demo scripts.
%



%% ==================== (BEGIN)

fprintf( '\n***** BEGIN (%s) *****\n\n', mfilename );


%% ==================== RUN DEMO SCRIPTS

% ---------- 2D demo

fprintf( ' > press any key to run the 2D demo\n' );
fprintf( '   ! figures and workspace variables will be cleared !\n' );
pause

run( 'demo_inversion_2d' );

% ---------- 3D demo (zero z-displacement)

fprintf( ' > press any key to run the 3D demo (zero z-displacement)\n' );
fprintf( '   ! figures and workspace variables will be cleared !\n' );
pause

run( 'demo_inversion_3d_z0' );

% ---------- 3D demo (sinusoidal z-displacement)

fprintf( ' > press any key to run the 3D demo (sinusoidal z-displacement)\n' );
fprintf( '   ! figures and workspace variables will be cleared !\n' );
pause

run( 'demo_inversion_3d_zsin' );


%% ==================== (END)

fprintf( '\n***** END (%s) *****\n\n', mfilename );



%%------------------------------------------------------------
%
% AUTHORS
%
%   Alexandros-Stavros Iliopoulos       ailiop@cs.duke.edu
%
% RELEASE
%
%   1.0.3 - March 21, 2019
%
% ------------------------------------------------------------
