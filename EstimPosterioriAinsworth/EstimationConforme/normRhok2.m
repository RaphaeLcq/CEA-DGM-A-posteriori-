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
% % Calcul de la norme de rho_K Article p1792 eq(45)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author : Raphael Lecoq CEA
%
%
% SYNOPSIS normRhoK2 = normRhok2(rho_K)
%
% GLOBALS:
%   - Aires(Nbtri,1)     : aire des triangles
%   - Nbtri              : nombre de triangle
% INPUT : 
%   - rho_K              : coefficients rho_n^(K)  dans l'article 1792 eq (42) 
% OUTPUT:
%   - normRhoK2(Nbtri,1) : norme L2 au carré de rho_K calculée avec la forme Article p1792 eq(45) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function normRhoK2 = normRhok2(rho_K)
  
global Nbtri 
global Aires 

%
loc_rho_K = zeros(3,2);
norm_rho_K2 = zeros(Nbtri,1);
L2contribution = zeros(Nbtri,1);
meanL2contribution = zeros(Nbtri,1);
%

for tri = 1:Nbtri 
  area48 = 48*Aires(tri); 
  loc_rho_K = squeeze(rho_K(tri,:,:)); %% rho_K(tri,:,:) : [1 x 3 x 2] --> squeeze --> [3 x 2]
  L2contribution(tri) = sum((loc_rho_K(1,:) + loc_rho_K(2,:)  + loc_rho_K(3,:)).^2)/area48;
  meanL2contribution(tri) = (sum(loc_rho_K(1,:).^2) + sum(loc_rho_K(2,:).^2) + sum(loc_rho_K(3,:).^2) )/area48;
endfor
%

normRhoK2 = L2contribution + meanL2contribution;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%  DEBUG   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global DEBUG
if (DEBUG)    
  for tri = 1:Nbtri
    rho1 = rho_K(tri, 1, :)/(2*Aires(tri));
    rho2 = rho_K(tri, 2, :)/(2*Aires(tri));
    rho3 = rho_K(tri, 3, :)/(2*Aires(tri));
    
    % Calcul des normes au carré pour chaque sommet
    norm2_rho1 = sum(rho1.^2);
    norm2_rho2 = sum(rho2.^2);
    norm2_rho3 = sum(rho3.^2);
    
    % Termes diagonaux : 2*Σ|ρK^(i)|²
    diagonal_terms = 2 * (norm2_rho1 + norm2_rho2 + norm2_rho3);
    
    % Termes croisés : Σ(i≠j) ρK^(i)*ρK^(j)
    cross_terms = 2*(sum(rho1 .* rho2) + ...
                       sum(rho1 .* rho3) + ...
                       sum(rho2 .* rho3));
    
    % Formule exacte de la norme L2 au carré
    norm2_DEBUG = (Aires(tri) / 12) * (diagonal_terms + cross_terms);
    
    if (abs(norm2_DEBUG - normRhoK2(tri)) > 1e-10)
      diff = norm2_DEBUG - normRhoK2(tri);
      fprintf("Probleme normRhok2 au triangle %i, la différence étant de %i\n",tri, diff);
      normRhoK2(tri) = norm2_DEBUG;
    endif
  endfor
  
endif

  

 
endfunction