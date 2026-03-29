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
% comp_estim1.m:
% Estimateur d'erreur issu de l'article
% "A POSTERIORI ERROR ESTIMATION FOR DISCONTINUOUS
% GALERKIN FINITE ELEMENT APPROXIMATION"
% MARK AINSWORTH
%
% SYNOPSIS [Estim_err, Estim_err_C, Estim_err_NC,Estim_err_data] = comp_estim1(Uh,idof)
%
% GLOBALS:
%   - CoorNeu(Nbpt,2)   : coordonnées (x, y) des sommets (noeuds P1)
%   - NumTri(Nbtri,3)   : liste des triangles (3 numéros de sommets)
%   - Nbtri             : nombre de triangles
%   - invLgEdg(Nbedg,1) : inverse de la longueur d'une face
%   - EdgNorm(Nbedg,:)  : vecteur normal non normalisé, orienté selon le maillage
%   - EdgUnit(Nbedg,:)  : vecteur normal normalisé, orienté selon le maillage
% INPUT :
%   - Uh                 : solution numérique
%   - idof(Nbtri,3)      : index des degrés de liberté
%
% OUTPUT:
%   - Estim_err                : estimateur d'erreur total
%   - Estim_err_C              : estimation d'erreur conforme
%   - Estim_err_NC             : estimation d'erreur de non conformité
%   - Estim_err_data           : estimation d'erreur d'oscillation des données
%   - Estim_err_C_without_curl : estimation d'erreur conforme sans l'amélioration avec le rotationnel eq (46) p1793
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Estim_err, Estim_err_C, Estim_err_NC,Estim_err_data,Estim_err_without_curl] = comp_estim1(Uh,idof)

global invLgEdg EdgUnit EdgNorm
global CoorNeu Nbtri NumTri
%%%%%%
global DEBUG
global visu

% Normalized vector EdgNorm
EdgUnit=EdgNorm.*invLgEdg;


%%%%%%%%%%%%%%%%%%% Estimateur non conforme défini article page 1795 eq (56)
[Estim_err_NC, Estim_err_NC2] = computeEstim_NC2(idof, Uh);

%%%%%%%%%%%%%%%%%%% Estimateur conforme
[Estim_err_C, Estim_err_C2, curl_rhok_norm2, norm_rho_K2, contribution_rhok2, norm2_Data, norm2_Neu] = compEstimCF2(idof,Uh);

%%%%%%%%%%%%%%%%%%% Estimation de l'oscillation des données
Estim_err_data2 = (sqrt(norm2_Data) + sqrt(norm2_Neu)).^2;
Estim_err_data = sqrt(sum(Estim_err_data2));

%%%%%%%%%%%%%%%%%%% Estimation sans correction par fonction bulle, pour comparer les résultats
Estim_err_C2_without_curl = (sqrt(norm_rho_K2) + sqrt(norm2_Data)  + sqrt(norm2_Neu)).^2 ;

%%%%%%%%%%%%%%%%%%% Estimateur des erreurs locales au carré
Estim_err_tot2 = Estim_err_C2 + Estim_err_NC2;
Estim_err_without_curl2 = Estim_err_C2_without_curl + Estim_err_C2;

%%%%%%%%%%%%%%%%%%% Estimateur d’erreur total
Estim_err = sqrt(sum(Estim_err_tot2));
Estim_err_without_curl = sqrt(sum(Estim_err_without_curl2));

%%%%%%%%%%%%%%%%%%% Erreur réelle
BrokenNorm2Real = GradL2Error(Uh, idof);

fprintf("Norme L2 de l'erreur non conforme = %.6e\n", Estim_err_NC);
fprintf("Norme L2  de la contribution curl rho_k = %.6e\n", sqrt(sum(curl_rhok_norm2)));
fprintf("Norme L2  de la contribution rho_k = %.6e\n", sqrt(sum(norm_rho_K2)));
fprintf("Norme L2 totale de rho_k = %.6e\n", sqrt(sum(contribution_rhok2)));
fprintf("Norme L2  de l'erreur conforme = %.6e\n", sqrt(sum(Estim_err_C2)));
fprintf("Norme L2  de l'oscillation des donnees = %.6e\n", Estim_err_data);

visualisationEstim(visu, Estim_err_NC2, curl_rhok_norm2, norm_rho_K2, contribution_rhok2, Estim_err_data2, Estim_err_C2_without_curl, Estim_err_C2, Estim_err_tot2, BrokenNorm2Real)
end


