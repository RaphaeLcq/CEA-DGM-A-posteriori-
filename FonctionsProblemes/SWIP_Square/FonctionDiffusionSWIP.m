%{
/****************************************************************************
* Copyright (c) 2022, CEA
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
% Author : Raphael Lecoq CEA
%
%
% SYNOPSIS %% %% Fonction f(r,theta) =  r^alpha ( a cos(alpha theta) + b*cos(alpha theta) ) avec a,b constants par morceaux
%
% Global - alpha 
%        - SWIP_coeffs(8,1) : coefficients (c_1,s_1, ..., c_4,s_4) pour les régions 1 à 4 dans le cas SWIP
%
% Inputs - rp(np,1) : Coordonnées polaire r 
%        - thetap(np,1) : Coordonnées polaire theta
% Output - fvals(np,1) : valeur de la fonction en ce point


function val = FonctionDiffusionSWIP(r,theta)
  global alpha;
  
  global SWIP_coeffs
  
  val = zeros(size(r));
  % Zone 1 : 0 < theta <= pi/2
  mask1 = (theta >= 0 & theta <= pi/2);
  val(mask1) =  r(mask1).^alpha .* ( SWIP_coeffs(1,1)*cos(alpha*theta(mask1)) + SWIP_coeffs(2,1) * sin(alpha*theta(mask1)));

  % Zone 2 : pi/2 < theta <= pi
  mask2 = (theta > pi/2 & theta < pi);
  val(mask2) =  r(mask2).^alpha .* ( SWIP_coeffs(3,1)*cos(alpha*theta(mask2)) + SWIP_coeffs(4,1)* sin(alpha*theta(mask2)));

  % Zone 3 : pi < theta <= 3*pi/2
  mask3 = (theta >= pi & theta <= 3*pi/2);
  val(mask3) =  r(mask3).^alpha .* ( SWIP_coeffs(5,1)*cos(alpha*theta(mask3)) + SWIP_coeffs(6,1)* sin(alpha*theta(mask3)));

  % Zone 4 : 3*pi/2 < theta <= 2*pi
  mask4 = (theta > 3*pi/2 & theta < 2*pi);
  val(mask4) =  r(mask4).^alpha .* ( SWIP_coeffs(7,1)*cos(alpha*theta(mask4)) + SWIP_coeffs(8,1)* sin(alpha*theta(mask4)));

  end
