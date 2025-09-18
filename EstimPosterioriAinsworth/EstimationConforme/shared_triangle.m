
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
% %% Given an edge of a triangle returns the quantities of interest on the edge of the triangle that shares the edge

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author : Raphael Lecoq CEA
%
%
% SYNOPSIS [shared_tri, shared_EdgeValue, shared_sigma, shared_muk]  = shared_triangle(Uh, idof, tri, xyNeu, edgeAGLO)
%
% GLOBALS:
%   - CoorBary(Nbtri,2)  : coordonnées des barycentres
%   - TriEdg(Nbtri,3)    : indices globaux des faces du triangle
%   - EdgTri(Nbedg,2)    : indices globaux des triangles qui partagent une face
%   - RefEdg(Nbedg,1)    : référence de Dirichlet/Neumann de la face 
%   - mu(Nbtri,3)        : mu(tri,i) indique si la normale de face i du triangle est dirigée vers l'extérieur du triangle (== 1) ou l'intérieur (== -1)
%   - kappa(Nbtri,1) : diffusion dans le triangle
% INPUT : 
%   - tri                : indice du triangle 
%   - Uh                 : solution numérique 
%   - idof(Nbtri,3)      : index des degrés de liberté
%   - edgeAGLO           : index global de la face partagée par tri et shared_tri
%   - xyNeu              : coordonnées des sommets de la face
% 
% OUTPUT:
%   - shared_tri         : triangle partagé par tri à la face de numérotation globale edgeAGLO
%   - shared_EdgeValue   : valeur de Uh à la face partagée dans le triangle partagé
%   - shared_sigma       :  valeur de sigma = kappa * grad(Uh) à la face partagée dans le triangle partagé
%   - shared_muk         : valeur de muk dans le triangle partagé
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [shared_tri, shared_EdgeValue, shared_sigma, shared_muk]  = shared_triangle(Uh, idof, tri, xyNeu, edgeAGLO, vertIGLO)

global EdgTri TriEdg RefEdg 
global mu kappa 
global DonneeDiri_P1

shared_tri = 0;
shared_EdgeValue = 0;
shared_sigma = zeros(2,1);
shared_muk = 0;


if (RefEdg(edgeAGLO) == 0) %% ddl interieur
  for i=1:2
    if (EdgTri(edgeAGLO,i) != tri)
      shared_tri = EdgTri(edgeAGLO,i); %% Global index of shared triangle
    endif
  endfor
  %% Shared triangle values
  shared_debT = idof(shared_tri,1); 
  shared_finT = idof(shared_tri,2);
  shared_UhT = Uh(shared_debT:shared_finT); %% Coefs of Uh in the shared triangle
  
  shared_EdgeValue = eval_fct_base(xyNeu,shared_tri)' * shared_UhT ; %% Uh(xyNeu) seen from shared triangle
  shared_sigma = kappa(shared_tri)*eval_grad_fct_base(xyNeu,shared_tri)' * shared_UhT ; %% sigma(xyNeu) seen from shared triangle
  shared_loc = find(TriEdg(shared_tri, :) == edgeAGLO); %find which local index is associated with the global index
  shared_muk = mu(shared_tri,shared_loc); %% mu(shared_triangle,edge)
  

elseif (RefEdg(edgeAGLO) > 0 & RefEdg(edgeAGLO) < 10) %% ddl Dirichlet
  shared_EdgeValue = DonneeDiri_P1(vertIGLO); %%L'estimateur utilise des données P1 (même si le RHS ne l'est pas);
  shared_muk = -1; 
endif


endfunction
