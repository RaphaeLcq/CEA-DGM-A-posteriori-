
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
%
%% Solution du problème  u(x,y) = x^3 - 3 x y^2
% Calcul des || \grad_f - \grad Uh ||_{ L^2(T) }^2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author :Raphael Lecoq CEA
%
%
% SYNOPSIS Calcule l'erreur du flux brisé pour le problème - Delta u = 0, u = x^3 - 3 x y^2 au bord dans le carré 
%
% GLOBALS:
%   - NumTri(Nbtri,3)        : liste des triangles (3 numéros de sommets)
%   - Nbtri                  : nombre de triangles
%   - CoorNeu(Nbpt,2)        : coordonnées des sommets
%   - Aires(Nbtri,1)         : aire des triangles
%   - npi                    : lambda*pi
%   - invDiaTri(Nbtri,1)     : inverse du diamètre du triangle
% INPUT : 
%   - Uh                     : solution numérique 
%   - idof(Nbtri,3)          : index des degrés de liberté
% 
% OUTPUT:
%   - NormGradDiff2(Nbtri,1) : norme L2 au carré de Grad( U_ex - U_h ) sur chaque triangle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function NormGradDiff2 = GradL2ErrorSquareNeumann(Uh, idof)
  global CoorNeu Nbtri NumTri invDiaTri Aires
  global npi 

  ndof = idof(Nbtri,2);
  NormGradDiff2 = zeros(Nbtri,1);
  NormDiff2 = zeros(Nbtri,1); %% CALCULER ERREUR L2


  for tri = 1:Nbtri
    IGLO = NumTri(tri,:); 
    CoorNeuT = CoorNeu(IGLO,:);
    aire = Aires(tri);
    invhT = invDiaTri(tri);
    deb = idof(tri,1); fin = idof(tri,2); ndofT = idof(tri,3);

    % Intégration sur triangle
    [xyp, wp, lambda, np] = IntTri_Ham7(CoorNeuT);
    awp = aire * wp';
    x = xyp(:,1); y = xyp(:,2);
    % Coordonnées barycentriques centrées
    grad = eval_grad_fct_base(xyp, tri);
    dPsiDGdx = squeeze( grad(:,1,:) ); %% transform ndof x 1 x np -> ndof x np 
    dPsiDGdy = squeeze( grad(:,2,:) );

    % Gradient Uh
    coeffs = Uh(deb:fin);  % coefficients locaux ndof x 1 
    gradUhX = coeffs' * dPsiDGdx;  % taille : 1 x np
    gradUhY = coeffs' * dPsiDGdy;


    % Somme partielle
    fx = 3*x.^2 - 3*y.^2;         
    fy = -6*x.*y;


    % Erreur quadratique
    diffX = fx' - gradUhX;
    diffY = fy' - gradUhY;
    diff2 = diffX.^2 + diffY.^2;

    NormGradDiff2(tri) = sum(awp .* diff2');
  endfor
  
endfunction