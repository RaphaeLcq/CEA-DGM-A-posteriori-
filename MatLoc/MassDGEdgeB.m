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
*
*****************************************************************************/
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author : Erell Jamelot CEA
%
% MassDGEdgeB.m:
% 
% DG-FEM (SIPG monomes P1 ou P2)
% 
% Matrice de stabilite pour une arete dans un triangle
%
% $\phi_i$ base locale dans $T_1$ (test fonction)
% $\psi_j$ base locale dans $T_2$ (trial function)
%
% Mxy(i,j)  = \int_F \phi_i \phi_j / |F|
% 
% SYNOPSIS Mxy=MassDGEdgeB(invhT,oT,CoorBaryT,CoorEdgF,CoorMilF)
%          
% INPUT : - invhT  : mise a l'echelle locale (inverse du diametre du triangle)
%         - oT : ordre d'approximation local au triangle
%         - CoorBaryT(1,2): coordonnees du barycentre du triangle
%         - CoorEdgF(2,2) : coordonnees de l'arete
%         - CoorMilF(1,2) : coordonnees du milieu de l'arete
% OUTPUT: - Mxy  : matrice de stabilite face
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Mxy=MassDGEdgeB(invhT,oT,CoorBaryT,CoorEdgF,CoorMilF)
%
% Schemas d'integration
[wp3,Mat3T,wp5,Mat5T]=ValPolyEdgeB(invhT,oT,CoorBaryT,CoorEdgF);
% Matrice de masse de l'arete, divisee par la longueur de l'arete
Mxy=MassDGEdge(invhT,oT,CoorBaryT,CoorMilF,wp3,Mat3T,wp5,Mat5T);