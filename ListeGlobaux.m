%{
/****************************************************************************
* Copyright (c) 2025, CEA
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

% GLOBAL - Nbtri  : global number of triangles
%        - Nbpt : global number of vertices
%        - Nbedg :  global number of edges

%        - CoorNeu(Nbpt,2) : coordinates (x, y) of vertices (nodes P1)
%        - CoorNeu2(Nbpt+Nbedg,2) : coordinates (x, y) of nodes P2
%        - CoorBary(Nbtri,2)   : coordinates of barycentres of triangles
%        - CoorMil(Nbedg,2)   : coordinates of the edges center

%        - LgEdg(Nbedg,1) : length of edges
%        - LgEdg2(Nbedg,1) : length of edges squared
%        - invLgEdg(Nbedg,1)  : inverse of the edge length
%        - DiaTri(Nbtri,1)     : diameter of triangles
%        - invDiaTri(Nbtri,2)   : inverse of diametres of triangles
%        - Aires(Nbtri,1) : area of triangles
%        - EdgNorm(Nbedg,2) : vector normal-edge, oriented tri1->tri2
%        - EdgUnit(Nbedg,2) : unitary vector normal-edge, oriented tri1->tri2


%        - NumTri(Nbtri,3) : list of triangles  (3 vertices number)
%        - NumTri2(4*Nbtri,3) : list of triangles of the mesh P2 (3 number of vertices)
%        - NumEdg(NbEdg,2) : number of 2 nodes of each vertice
%        - NumEdgB(NbEdgB,2) : number of 2 nodes of each vertice
%		     - TriEdg(Nbtri,3) : for each triangle, TriEdg(l,i) is the number of vertice opposite from vertice NumTri(l,i)
%                  (3 number of edges - full matrix Nbtri x 3)
%		     - EdgTri(Nbedg,2) : for each vertice, EdgTri(a,:) gives the number of the 2 triangles sharing vertice
%                                 EdgTri(a,2) = 0 if a is on the domain border
%		     - SomOpp(NbEdg,2) : number of the opposite vertice to the edge in each sharing triangle
%                                  SomOpp(a,2) = 0 if a is on the domain border

%        - RefTri(Nbtri,1) : reference of triangles of the mesh
%        - RefTri2(Nbtri2,1) : reference of triangles of the mesh raffine
%        - RefNeu(Nbpt,1) : reference of vertices
%		     - RefEdg(Nbedg,1) : reference of each vertice
%		     - RefEdgB(NbEdgB,1) : reference of each vertice (boundaries)
%        - RefNeu2(Nbpt2,1) : reference of vertices followed of the middle of the edges

%        - kappa(Nbtri,1)   : diffusion in the triangle
%        - mu(Nbtri,3)        : mu(tri,i) indicate if the normal of face i of the triangle is directed tower the outside of the triangle (== 1) or its interior (== -1)
%        - sigTri(Nbtri,1)    : sigma of the triangle
%        - etaEdg(NbEdg,1) : value of eta on each edge

%        - ordre  :  approximation order
%        - fig : number of figure if visualisation
%        - mi  : number of the mesh
%        - npi : lambda * pi

%        - visu : chosen visualisation by the user in mainPoissonDG

