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
% Author : Raphael Lecoq CEA
%
%
% SYNOPSIS %%%%%%% Projection d'une fonction lisse u_D sur les fonctions P1 par morceaux CONTINUES sur le bord uD_P1

%

% GLOBAL - Nbtri  : nombre de triangle
%		     - EdgTri(Nbedg,2) : Pour chaque arete, EdgTri(a,:) donne les numeros des 2 triangles de chaque arete 
%                                 EdgTri(a,2) = 0 si a est sur la frontiere
%		     - TriEdg(Nbtri,3) : Pour chaque triangle, TriEdg(l,i) est le numero de l'arete opposee au sommet NumTri(l,i)
%                  (3 numeros des aretes - matrice entiere Nbtri x 3)
%        - NumTri(Nbtri,3) : liste de triangles  (3 numeros de sommets) 
%        - NumEdg(NbEdg,2) : Numero des 2 noeuds de chaque arete
%        - CoorNeu(Nbpt,2) : coordonnees (x, y) des sommets (noeuds P1)
%        - Nbpt : nombre de sommets
%        - Nbedg : nombre d'aretes
%		     - RefEdg(Nbedg,1) : Reference de chaque arete 

%

function projP1c_Edge()
  
global EdgTri TriEdg NumTri NumEdg CoorNeu Nbpt Nbedg RefEdg
global LgEdg
global DonneeDiri_P1
% Initialiser DonneeDiri_P1 avec des zéros pour tous les sommets
DonneeDiri_P1 = sparse(Nbpt, 1);

%% Sommets globaux associés à un edge Dirichlet global
IndDirichlet = find(RefEdg > 0 & RefEdg < 10);
NbDiri = length(IndDirichlet);

% Si aucun bord de Dirichlet, retourner avec DonneeDiri_P1 = 0 partout
if NbDiri == 0
    return;
end

% Créer la liste des sommets Dirichlet uniques et leur mapping
sommets_diri = unique(NumEdg(IndDirichlet, :));
NbSommets = length(sommets_diri);

% Créer un mapping : sommet global -> indice local dans le système
sommet_to_local = sparse(max(sommets_diri), 1);
for k = 1:NbSommets
    sommet_to_local(sommets_diri(k)) = k;
end

% Matrices et vecteur du système local
M = sparse(NbSommets, NbSommets);
b = sparse(NbSommets, 1);

for edge = IndDirichlet'
  i_global = NumEdg(edge,1);
  j_global = NumEdg(edge,2);
  
  % Conversion vers indices locaux
  i_local = sommet_to_local(i_global);
  j_local = sommet_to_local(j_global);
  
  L = LgEdg(edge); % longueur de l'arête
  % Masse locale
  Medg = L/6 .* [2 1; 1 2];
  
  % Assemblage avec indices locaux
  M(i_local,i_local) = M(i_local,i_local) + Medg(1,1);
  M(i_local,j_local) = M(i_local,j_local) + Medg(1,2);
  M(j_local,i_local) = M(j_local,i_local) + Medg(2,1);
  M(j_local,j_local) = M(j_local,j_local) + Medg(2,2);
  
  XY = [CoorNeu(i_global,:); CoorNeu(j_global,:)];
  [xyp, wp, lambda, np] = IntEdg_Boo5(XY);
  uD_quad = eval_donneesDirichlet(xyp(:,1),xyp(:,2));
  awp = LgEdg(edge) .* wp';
  
  % Intégrale locale RHS avec indices locaux
  b(i_local) += sum(awp .* lambda(:,1) .* uD_quad);
  b(j_local) += sum(awp .* lambda(:,2) .* uD_quad);
end

% Résolution du système local
uD_P1_local = M \ b;

% Reconstruction du vecteur global : 
% DonneeDiri_P1(s) = 0 si s n'est pas sur un bord de Dirichlet
% DonneeDiri_P1(s) = valeur calculée si s est sur un bord de Dirichlet
for k = 1:NbSommets
    s = sommets_diri(k);  % indice global du sommet
    DonneeDiri_P1(s) = uD_P1_local(k);  % Assigner la valeur calculée
end

endfunction