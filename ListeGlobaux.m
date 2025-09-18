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
%
% SYNOPSIS
%

% GLOBAL - Nbtri  : nombre de triangle
%        - Nbpt : nombre de sommets
%        - Nbedg : nombre d'aretes

%        - CoorNeu(Nbpt,2) : coordonnees (x, y) des sommets (noeuds P1)
%        - CoorNeu2(Nbpt+Nbedg,2) : coordonnees (x, y) des noeuds P2
%        - CoorBary(Nbtri,2)   : Coordonnees des barycentres de triangles
%        - CoorMil(Nbedg,2)   : Coordonnees des milieux d'aretes

%        - LgEdg(Nbedg,1) : longueurs des aretes
%        - LgEdg2(Nbedg,1) : longueurs des aretes au carre
%        - invLgEdg(Nbedg,1)  : inverse de la longueur de la face
%        - DiaTri(Nbtri,1)     : diamètre du triangle
%        - invDiaTri(Nbtri,2)   : inverse des diametres des triangles
%        - Aires(Nbtri,1) : aires des triangles
%        - EdgNorm(Nbedg,2) : vecteurs face-normale, orientes tri1->tri2
%        - EdgUnit(Nbedg,2) : vecteurs face-normale unitaires, orientes tri1->tri2


%        - NumTri(Nbtri,3) : liste de triangles  (3 numeros de sommets) 
%        - NumTri2(4*Nbtri,3) : liste de triangles du maillage P2 (3 numeros de sommets)
%        - NumEdg(NbEdg,2) : Numero des 2 noeuds de chaque arete
%        - NumEdgB(NbEdgB,2) : Numero des 2 noeuds de chaque arete
%		     - TriEdg(Nbtri,3) : Pour chaque triangle, TriEdg(l,i) est le numero de l'arete opposee au sommet NumTri(l,i)
%                  (3 numeros des aretes - matrice entiere Nbtri x 3)
%		     - EdgTri(Nbedg,2) : Pour chaque arete, EdgTri(a,:) donne les numeros des 2 triangles de chaque arete 
%                                 EdgTri(a,2) = 0 si a est sur la frontiere
%		     - SomOpp(NbEdg,2) : Numero du sommet oppose a l'arete dans chaque triangle
%                                  SomOpp(a,2) = 0 si a est sur la frontiere 

%        - RefTri(Nbtri,1) : reference des triangles du maillage 
%        - RefTri2(Nbtri2,1) : reference des triangles du maillage raffine
%        - RefNeu(Nbpt,1) : reference des sommets
%		     - RefEdg(Nbedg,1) : Reference de chaque arete 
%		     - RefEdgB(NbEdgB,1) : Reference de chaque arete
%        - RefNeu2(Nbpt2,1) : reference des sommets suivis des milieux d'aretes

%        - kappa(Nbtri,1)   : diffusion dans le triangle 
%        - mu(Nbtri,3)        : mu(tri,i) indique si la normale de face i du triangle est dirigée vers l'extérieur du triangle (== 1) ou l'intérieur (== -1)
%        - sigTri(Nbtri,1)    : sigma du triangle (je sais pas ce que c'est, on l'utilise dans EtaParam.m, peut être le nombre de faces ?) 

%        - ordre  : ordre d'approximation
%        - fig : numero de figure si visualisation
%        - mi  : numero du maillage
%        - npi : lambda * pi

