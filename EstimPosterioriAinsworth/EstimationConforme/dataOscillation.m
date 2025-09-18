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
% Author :Raphael Lecoq  CEA
%
% Calcul des quantités faisant intervenir rho_K dans l'article p1794 eq (53)
%
% SYNOPSIS [norm2_Data, norm2_Diri, norm2_Neu] = dataOscillation()
%
% GLOBALS:
%   - DiaTri(Nbtri,1)     : diamètre du triangle
%   - Nbtri               : nombre de triangle
%   - NumTri(Nbtri,3)     : indices globaux des sommets du triangle
%   - CoorNeu(Nbtri,3)    : coordonnées des sommets du triangle
%   - Aires(Nbtri,1)      : aire d'un triangle
%   - kappa(Nbtri,1)         : diffusion dans le triangle 
% INPUT : 
% 
% OUTPUT:
%   - norm2_Data(Nbtri,1) : carré de la norme L2 |f - average(f)_T|
%   - norm2_Diri(Nbtri,1) : carré de la norme L2 |g_D - average(g_D)_F| où g_D est la valeur de Uh aux bords de Dirichlet 
%   - norm2_Neu(Nbtri,1)  : carré de la norme L2 |g_N - average(g_N)_F| où g_N est la valeur de Grad(Uh) aux bords de Neumann 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [norm2_Data, norm2_Neu] = dataOscillation()

global DiaTri Aires LgEdg
global NumTri Nbtri Nbpt NumEdg CoorNeu kappa RefEdg EdgTri

norm2_Data = zeros(Nbtri,1);
norm2_Neu = zeros(Nbtri,1);




Cp = DiaTri./pi; %% Constante de Poincaré Theorem 4 eq (61) p1797

%%%%% Calcul de l'oscillation du terme source f(x,y) 
for tri = 1:Nbtri
  area = Aires(tri);
  XY = CoorNeu(NumTri(tri,:), :);
  [xyp, wp, ~, ~] = IntTri_Ham7(XY);
  awp = area * wp';  
  %
  fvals = evalDonnees_Source(xyp(:,1),xyp(:,2));  
  mean_f = sum(awp .* fvals) / area;         % Moyenne de f sur le triangle
  diff = fvals - mean_f;
  norm2_Data(tri) = sum(awp .* diff.^2);     % Intégrale de (f - moyenne)^2 sur le triangle

endfor 
norm2_Data = Cp.^2.*norm2_Data./kappa;  


%%%%%%% Calcul de l'oscillation des données sur le bords de Neumann
IndNeu = find(RefEdg >= 10);
IndNeu = IndNeu';
EdgeValue = zeros(Nbpt,1);
for edge = IndNeu
  vertices = NumEdg(edge,:);
  XY = CoorNeu(vertices,:);
  [xyp, wp, ~, ~] = IntEdg_Boo5(XY);
  awp = LgEdg(edge) * wp';
  
  PhiNeumann = eval_donneesNeumann(xyp(:,1), xyp(:,2),edge);
  mean_PhiNeumann = sum(awp.* PhiNeumann) / LgEdg(edge);
  
  diff_Diri = PhiNeumann - mean_PhiNeumann;
  
  Norm2_Edge = sum(awp.* diff_Diri.^2);
  tri = EdgTri(edge,1);
  [EdgeDiam,EdgeHeight] = GeoQuantities(tri,edge); %% Pour calculer Ct de trace Lemme 11 eq (64) p1798
  norm2_Neu(tri) = 2*Cp(tri)*(Cp(tri) + EdgeDiam)/EdgeHeight * Norm2_Edge / kappa(tri);
endfor
  
endfunction
