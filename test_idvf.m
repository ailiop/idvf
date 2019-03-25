% test_idvf.m
%
% Test idvf functionality by running the included demo scripts.
%
% ----------------------------------------------------------------------
%
% Copyright (C) 2019, Department of Computer Science, Duke University
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
%   1.0.3 - March 25, 2019
%
% ------------------------------------------------------------
