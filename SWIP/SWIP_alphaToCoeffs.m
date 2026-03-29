%{
/****************************************************************************
* Copyright (c) 2025, CEA
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
* OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author : Raphael Lecoq, CEA
%
%
% SYNOPSIS
% Compute the coefficient of diffusion D chosen for the square SWIP problem.
%
% GLOBAL : None
% INPUT :
%   - alpha        : power of the radial coordinate for the square SWIP problem
% OUTPUT:
%   - D    : value of the diffusion for the square SWIP problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [D] = SWIP_alphaToCoeffs(alpha)

  x_1 = cos(alpha*pi/2);
  y_1 = sin(alpha*pi/2);
  if (abs(y_1) < 1e-16) %% M is always non-invertible, which means any D != 1 will do
    D = 0.5;
    disp('yahou ');
  else
    z_1 = (x_1/y_1);
    z_1_pwr2 = (x_1/y_1)^2;
    D_plus = (2*z_1_pwr2 + 1) + 2*z_1*sqrt(z_1_pwr2 + 1);
    D_min = (2*z_1_pwr2 + 1) - 2*z_1*sqrt(z_1_pwr2 + 1);
    D = min(D_min, D_plus);
  endif
  M = zeros(8,8);
  A = zeros(5,2,2);
  B = zeros(5,2,2);
  for n = 0:4
    theta = n*pi/2;
    x_n = cos(alpha*theta);
    y_n = sin(alpha*theta);
    A(n+1,:,:) = [ x_n , y_n ; -y_n , x_n ];
    B(n+1,:,:) = -[ x_n , y_n ; -D*y_n , D*x_n ];
  endfor
  M(1:2,1:2) = squeeze(A(2,:,:));
  M(1:2,3:4) = squeeze(B(2,:,:));
  M(3:4,3:4) = squeeze(B(3,:,:));
  M(3:4,5:6) = squeeze(A(3,:,:));
  M(5:6,5:6) = squeeze(A(4,:,:));
  M(5:6,7:8) = squeeze(B(4,:,:));
  M(7:8,1:2) = squeeze(A(1,:,:));
  M(7:8,7:8) = squeeze(B(5,:,:));
  global SWIP_coeffs =  robust_null_combined(M)

  SWIP_coeffs = SWIP_coeffs(:,1);
  global DiffusionConstant = D
  fprintf('Diffusion for alpha = %i must be D = %i\n',alpha,D)
endfunction
