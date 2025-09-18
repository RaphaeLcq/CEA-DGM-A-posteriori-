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
% Author : Raphaël Lecoq CEA
%
%
% SYNOPSIS Affiche la valeur des moyennes pondérées au niveau des faces pour la méthode SWIP
%

% GLOBAL - NumEdg(NbEdg,2) : Numero des 2 noeuds de chaque arete
%        - CoorNeu(Nbpt,2) : coordonnees (x, y) des sommets (noeuds P1)
%		     - RefEdg(Nbedg,1) : Reference de chaque arete 
%        - Nbedg : nombre d'aretes
%		     - EdgTri(Nbedg,2) : Pour chaque arete, EdgTri(a,:) donne les numeros des 2 triangles de chaque arete 
%                                 EdgTri(a,2) = 0 si a est sur la frontiere
%        - kappa(Nbtri,1)   : diffusion dans le triangle 


function plotAvgKappa()
    global NumEdg CoorNeu RefEdg Nbedg EdgTri kappa

    if isempty(NumEdg) || isempty(CoorNeu) || isempty(RefEdg) || isempty(Nbedg) || isempty(kappa)
        error('Variables globales manquantes');
    end

    figure;
    hold on;

    for edge = 1:Nbedg
        % Coordonnées des sommets de l'arête
        s1 = NumEdg(edge, 1);
        s2 = NumEdg(edge, 2);
        x_coords = [CoorNeu(s1, 1), CoorNeu(s2, 1)];
        y_coords = [CoorNeu(s1, 2), CoorNeu(s2, 2)];

        % Tracer l'arête en noir
        plot(x_coords, y_coords, 'k-', 'LineWidth', 1);

        % Si arête interne
        if RefEdg(edge) == 0 
            tri1 = EdgTri(edge,1);
            tri2 = EdgTri(edge,2);
            [AvgK1,AvgK2] = HarmonicAverageKappa(edge);
            kappa1 = kappa(tri1);
            kappa2 = kappa(tri2);
            fprintf('Tri1 = %i et Tri2 = %i, kappa1 = %i et kappa2 = %i, AvgK = [%i,%i]\n',tri1,tri2,kappa1,kappa2,AvgK1,AvgK2)
            x_mid = mean(x_coords);
            y_mid = mean(y_coords);
            text(x_mid, y_mid, sprintf('%.3f', AvgK1), ...
                'Color', 'r', 'FontSize', 8, 'HorizontalAlignment', 'center');
        end
    end

    title('Moyenne harmonique de \kappa uniquement aux interfaces de contraste');
    axis equal;
    hold off;
end
