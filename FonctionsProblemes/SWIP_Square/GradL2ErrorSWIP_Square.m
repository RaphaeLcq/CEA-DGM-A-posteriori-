
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
%% Solution du problème  f(r,theta) = r^(alpha)( a cos(theta alpha) + b sin(theta alpha) )
% Calcul des || \grad_f - \grad Uh ||_{ L^2(T) }^2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author : Raphael Lecoq CEA
%
%
% SYNOPSIS Calcule l'erreur du flux brisé pour SWIP dans un carré partitionné en 4 carrés
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
%   - NormGradDiff2(Nbtri,1) : norme L2 au carré de Grad(kappa^1/2 * U_ex - U_h ) sur chaque triangle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function NormGradDiff2 = GradL2ErrorSWIP_Square(Uh, idof)
  global CoorNeu Nbtri NumTri invDiaTri Aires
  global RefTri
  global SWIP_coeffs alpha 

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


    region = RefTri(tri);
    [rp,thetap] = CartesianToPolarCentered(xyp(:,1),xyp(:,2));
      % Préallocation
    gradx = zeros(size(rp));
    grady = zeros(size(rp));

    % --- Fonction de base
    % Chaque zone a une combinaison linéaire différente
    % Définir coefficients selon zone
    coeff_cos = zeros(size(rp));
    coeff_sin = zeros(size(rp));

    mask1 = (thetap >= 0 & thetap <= pi/2);
    coeff_cos(mask1) = SWIP_coeffs(1,1);
    coeff_sin(mask1) = SWIP_coeffs(2,1);

    mask2 = (thetap > pi/2 & thetap < pi);
    coeff_cos(mask2) = SWIP_coeffs(3,1);
    coeff_sin(mask2) = SWIP_coeffs(4,1);

    mask3 = (thetap >= pi & thetap <= 3*pi/2);
    coeff_cos(mask3) = SWIP_coeffs(5,1);
    coeff_sin(mask3) = SWIP_coeffs(6,1);

    mask4 = (thetap > 3*pi/2 & thetap < 2*pi);
    coeff_cos(mask4) = SWIP_coeffs(7,1);
    coeff_sin(mask4) = SWIP_coeffs(8,1);

    % f(r,theta) = r^alpha * ( coeff_cos * cos(alpha*theta) + coeff_sin * sin(alpha*theta) )
    base = coeff_cos .* cos(alpha*thetap) + coeff_sin .* sin(alpha*thetap);

    % --- Dérivées polaires
    df_dr = alpha * rp.^(alpha-1) .* base;
    df_dtheta = rp.^alpha .* ( -alpha*coeff_cos .* sin(alpha*thetap) + alpha*coeff_sin .* cos(alpha*thetap) );

    % --- Passage en cartésien
    fx = cos(thetap) .* df_dr - (sin(thetap)./rp) .* df_dtheta;
    fy = sin(thetap) .* df_dr + (cos(thetap)./rp) .* df_dtheta;



    % Erreur quadratique
    diffX = fx' - gradUhX;
    diffY = fy' - gradUhY;
    diff2 = diffX.^2 + diffY.^2;

    NormGradDiff2(tri) = sum(awp .* diff2');
  endfor
  
endfunction