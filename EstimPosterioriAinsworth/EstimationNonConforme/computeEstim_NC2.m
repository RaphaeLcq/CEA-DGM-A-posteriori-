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
% Author : CEA
%
% Calcul de l'erreur de non conformité page 1795 eq (56)
% 
% SYNOPSIS [Estim_err_NC, Estim_err_NC2]= computeEstim_NC2(idof, Uh)
%
% GLOBALS:
%   - Nbtri                  : Nombre de triangles
%   - CoorBary(Nbtri,2)      : coordonnées (x, y) des barycentres (noeuds P0)
%   - CoorMil(Nbedg,2)       : coordonnées (x, y) des milieux (noeuds P1 dG)
%   - TriEdg(Nbtri,2)        : indices globaux des faces du triangle
%   - Aires(Nbtri,1)         : aire du triangle
%   - invDiaTri(Nbtri,1)     : inverse de l'aire du triangle
%   - ordre                  : ordre de la méthode numérique (P1 ou P2)
% INPUT : 
%   - idof(Nbtri,3)          : index des degrés de liberté
%   - Uh                     : solution numérique 
% 
% OUTPUT:
%   - Estim_err_NC           : estimateur d'erreur non conforme
%   - Estim_err_NC2(Nbtri,1) : carré de l'estimateur d'erreur local pour chaque triangle 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Estim_err_NC, Estim_err_NC2]= computeEstim_NC2(idof, Uh)
  
global Nbtri CoorBary CoorMil TriEdg
global Aires invDiaTri 
global ordre diffusion

diffusion = ones(Nbtri,1); %% heterogeneous diffusion .. or not 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U_star = projCoorBaryP1(Uh, idof); %% U_star = sum_n Un^* \lambda_n where \lambda_n is the n-th barycentric coordinate
                                   % Article p1795 eq (56)
grad_U_star = gradCoorBary(idof,U_star); %% Article p1796 eq (57)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ù


%% Norme |a * grad (u - u^*)|^2 Article p1796 eq (57)
Estim_err_NC2 = zeros(Nbtri,1); 
Estim_err_NC = 0;
%
aThird = 1/3; %% pondération pour la quadrature P2
n_quadrature = 2*ordre - 1; %% Nbre points de quadrature, 1 si ordre == 1, 3 si ordre == 2
%

for tri = 1:Nbtri 
  
    area = Aires(tri);
    grad_U_star_tri = grad_U_star(tri, :);
    debT = idof(tri,1); finT = idof(tri,2);
    Uh_T = Uh(debT:finT); 
    
    %%%%% Quadrature pour norme L2
    if ordre == 1
      xq = CoorBary(tri, :); % quadrature P0      
      awq = area ; 
    else %% ordre == 2
      xq = CoorMil(TriEdg(tri,:),:); % quadrature P2 
      awq = area * aThird; 
    endif
    
    for i = 1:n_quadrature
      grad_diff = eval_grad_fct_base(xq(i,:), tri)' * Uh_T - grad_U_star_tri' ;
      Estim_err_NC2(tri) += diffusion(tri) * awq * sum( grad_diff.^2, 1);
    endfor
end
%

Estim_err_NC = sqrt(sum(Estim_err_NC2)); %%Article p 1783

  endfunction