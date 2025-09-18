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
% Author : Raphaël Lecoq  CEA
%
%
% SYNOPSIS Affiche la référence des faces 
%

% GLOBAL - NumEdg(NbEdg,2) : Numero des 2 noeuds de chaque arete
%        - CoorNeu(Nbpt,2) : coordonnees (x, y) des sommets (noeuds P1)
%		     - RefEdg(Nbedg,1) : Reference de chaque arete 
%        - Nbedg : nombre d'aretes


function plotEdgeRef()
    global NumEdg CoorNeu RefEdg Nbedg

    if isempty(NumEdg) || isempty(CoorNeu) || isempty(RefEdg) || isempty(Nbedg)
        error('Variables globales manquantes');
    end

    figure;
    hold on;

    for e = 1:Nbedg
        % Sommets de l'arête
        s1 = NumEdg(e,1);
        s2 = NumEdg(e,2);
        X = [CoorNeu(s1,1), CoorNeu(s2,1)];
        Y = [CoorNeu(s1,2), CoorNeu(s2,2)];

        % Tracer l'arête
        plot(X, Y, 'k-');

        % Milieu
        xm = mean(X);
        ym = mean(Y);

        % Référence de l'arête
        ref = RefEdg(e);

        % Affichage
        text(xm, ym, num2str(ref), ...
             'Color', 'm', 'FontSize', 8, 'HorizontalAlignment', 'center');
    end

    title('Références des faces/arêtes');
    axis equal;
    hold off;
end
