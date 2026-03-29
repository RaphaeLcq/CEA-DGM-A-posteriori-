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
%%%%%%%%%%%%% Projection d'une fonction dG sur une fonction P1 LG par moyenne aux sommets
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author : CEA
%
% 
% SYNOPSIS projUh = projCoorBaryP1(Uh, idof)
%
% GLOBALS:
%   - Nbpt                   : nombre de sommets
%   - CoorNeu(Nbpt,2)        : coordonnées (x, y) des sommets (noeuds P1 dG)
%   - RefNeu(Nbpt,1)         : référence de Dirichlet du sommet
%   - Aires(Nbtri,1)         : aire du triangle
%   - invDiaTri(Nbtri,1)     : inverse du diamètre du triangle
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


% 
%
%
function projUh = projCoorBaryP1(Uh, idof)
  

global CoorNeu NumTri Nbpt RefNeu NumEdg RefEdg

projUh = zeros(Nbpt,1); %% liste des valeurs aux sommets
global DonneeDiri_P1; %% Function de Dirichlet aux sommets
    
for s = 1:Nbpt
  is_in_edge = (NumEdg(:,1) == s) | (NumEdg(:,2) == s);
  edges = find(is_in_edge);
  RefDiri = find(RefEdg(edges) > 0 & RefEdg(edges) < 10);
  VerticeOnDiri = 1-isempty(RefDiri);
  if (VerticeOnDiri==0)% ddl interieur ou Neumann
    Ind=find(sum(NumTri==s,2)); % triangles ayant pour sommet s
    stencil=length(Ind); % nombre de triangles ayant pour sommet s
    xk = CoorNeu(s,:); % coordonnées de s
    
    for tri = transpose(Ind)
      %% base des phi_i évaluée en xk
      val_fct_base = eval_fct_base(xk, tri); 
      % coefficients locaux dans Uh
      debT = idof(tri,1); 
      finT = idof(tri,2);
      Uh_T = Uh(debT:finT); %% coefs dans la base DG tq Uh = sum_i u_i phi_i de taille (n x 1)
      % évaluation
      projUh(s) += val_fct_base' * Uh_T; %% = sum_i u_i phi_i(x,y) = U_h(x,y)
      
    endfor
    projUh(s) = projUh(s) / stencil;
    
  elseif (VerticeOnDiri) %ddl dirichlet, projUh(s) = f_D(xk) si f_D est la fonction de dirichlet 
      projUh(s) = DonneeDiri_P1(s);
  endif
  
endfor

%%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global DEBUG visu
if (DEBUG && visu == 1)
  fig = 75;
  global NumTri CoorNeu
  NT=NumTri; CN=CoorNeu;
  figure(fig);
  trisurf(NT,CN(:,1),CN(:,2),projUh);
  view(2);
  colormap ("jet");
  shading interp
  title("Proj Bary P1")
  colorbar;
endif
endfunction 