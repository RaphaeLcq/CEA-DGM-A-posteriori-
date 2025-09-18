
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
%% Solution du problème harmonique dans le carré f(r,theta) = r^(alhpa) * (asin(alpha theta) + b cos(alpha theta))
% Calcul des || \grad_f - \grad Uh ||_{ L^2(T) }^2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author : Raphael Lecoq CEA
%
%
% SYNOPSIS Calcule l'erreur du flux brisé pour le problème harmonique dans le carré unitaire
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


function NormGradDiff2 = GradL2ErrorHarmonicSquare(Uh, idof)
  global CoorNeu Nbtri NumTri invDiaTri Aires
  global npi 

  ndof = idof(Nbtri,2);
  NormGradDiff2 = zeros(Nbtri,1);
  NormDiff2 = zeros(Nbtri,1); %% CALCULER ERREUR L2
  
  global HarmonicSquare_alpha;
  global HarmonicSquare_a ;
  global HarmonicSquare_b ;
  
  AlphaMinus1 = HarmonicSquare_alpha - 1;
  
  for tri = 1:Nbtri
    IGLO = NumTri(tri,:); 
    CoorNeuT = CoorNeu(IGLO,:);
    aire = Aires(tri);
    invhT = invDiaTri(tri);
    debT = idof(tri,1); finT = idof(tri,2); 

    % Intégration sur triangle
    [xyp, wp, lambda, np] = IntTri_Ham7(CoorNeuT);
    awp = aire * wp';

    % Coordonnées barycentriques centrées
    grad = eval_grad_fct_base(xyp, tri);
    dPsiDGdx = squeeze( grad(:,1,:) ); %% transform ndof x 1 x np -> ndof x np 
    dPsiDGdy = squeeze( grad(:,2,:) );

    % Gradient Uh
    Uh_T = Uh(debT:finT);  % coefficients locaux ndof x 1 
    gradUhX = Uh_T' * dPsiDGdx;  % taille : 1 x np
    gradUhY = Uh_T' * dPsiDGdy;

    % Gradient f
    [rp, thetap] = CartesianToPolarCentered(xyp(:,1), xyp(:,2));
    df_dr     = HarmonicSquare_alpha .* rp.^AlphaMinus1 .*( HarmonicSquare_a .* sin(AlphaMinus1.*thetap) + HarmonicSquare_b .* cos(AlphaMinus1 .* thetap) );
    df_dtheta = HarmonicSquare_alpha .* rp.^AlphaMinus1 .*( HarmonicSquare_a .* cos(AlphaMinus1.*thetap) - HarmonicSquare_b .* sin(AlphaMinus1 .* thetap) );

    
    % Erreur quadratique
    diffX = df_dr' - gradUhX;
    diffY = df_dtheta' - gradUhY;

    diff2 = diffX.^2 + diffY.^2;

    NormGradDiff2(tri) = sum(awp .* diff2');
  endfor
  
endfunction