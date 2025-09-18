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
%% Évaluation de la valeur du gradient des fonctions de base d'un triangle donné à des coordonnées globales données 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author : Andrew Peitavy, Raphael Lecoq CEA
%
% SYNOPSIS  [val] = eval_fct_base(xk,tri)
%
% GLOBALS:
%
% INPUT : 
%   - xk : points d'évaluation des fonctions de bases
%   - tri  : numérotation du triangle
% OUTPUT:
%   - val : valeurs du gradient des fonctions de bases évalués en xk
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [val] = eval_grad_fct_base(xk, tri)

    global ordre CoorBary invDiaTri
    CoorBary_k=CoorBary(tri,:);
    invDiaTri_k=invDiaTri(tri,:);
    
    size=size(xk);
    if(size(2)!=2)
      disp("Erreur, mauvaise dimension donné à la fonction eval_fct_base");
    endif 

    % Initialisation du tableau du gradient des fonctions de base
    if ordre == 1
        val = zeros(3,2,size(1));
    elseif ordre == 2
        val = zeros(6,2,size(1));
    end
    
    for nb_pts=1:size(1)  
    % Coordonnées locales par rapport au barycentre

    x_1 = (xk(nb_pts,1) - CoorBary_k(1))*invDiaTri_k;
    x_2 = (xk(nb_pts,2) - CoorBary_k(2))*invDiaTri_k;

    % Gradient des fonctions de base du premier ordre
    val(1, :,nb_pts) = [0, 0];
    val(2, :,nb_pts) = [invDiaTri_k, 0];
    val(3, :,nb_pts) = [0, invDiaTri_k];

    % Gradient des fonctions supplémentaires pour le second ordre
    if ordre == 2
        val(4, :,nb_pts) = [invDiaTri_k * x_2, invDiaTri_k * x_1];
        val(5, :,nb_pts) = [2 * invDiaTri_k * x_1, 0];
        val(6, :,nb_pts) = [0, 2 * invDiaTri_k * x_2];
    end
    
    endfor
end