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
% SYNOPSIS Affiche les triangles où l'équilibre des flux n'est pas respecté dans l'article d'Ainsworth 
%

% GLOBAL 
% INPUTS:
%   - idof : index des degrés de liberté
%   - Uh : solution numérique
%   - tolerance : seuil d'erreur (défaut: 1e-10)



function [erreur_triangles, erreur_valeurs] = VisualiseDEBUGTriangles(idof, Uh, tolerance)
% Visualise les triangles où le test DEBUG n'est pas vérifié

global NumTri TriEdg EdgTri CoorNeu CoorBary RefEdg
global etaEdg invDiaTri invLgEdg mu EdgUnit kappa
global Aires npi npi2 LgEdg Lshape
global Nbtri

if nargin < 3
    tolerance = 1e-10;
end

% Initialisation
erreur_triangles = [];
erreur_valeurs = [];
int_gk_values = zeros(Nbtri, 1);
int_f_values = zeros(Nbtri, 1);
erreur_equilibre = zeros(Nbtri, 1);
triangle_type = zeros(Nbtri, 1); % 0: intérieur, 1: bord Dirichlet, 2: autre bord

fprintf('=== VERIFICATION DEBUG POUR TOUS LES TRIANGLES ===\n');

% Boucle sur tous les triangles
for tri = 1:Nbtri
    % Calcul des fonctions de flux gk pour ce triangle
    gk = fluxFunctions(tri, idof, Uh); % Version silencieuse
    AGLO = TriEdg(tri,:);
    IGLO = NumTri(tri,:);
    CoorNeuT = CoorNeu(IGLO,:);
    
    % Calcul de l'intégrale sur le bord: int_{\partial K} gk
    int_gk_edge = 0;
    for i_edge_loc = 1:3
        int_gk_edge = int_gk_edge + LgEdg(AGLO(i_edge_loc)) * 0.5 * sum(gk(i_edge_loc,:));
    end
    int_gk_values(tri) = int_gk_edge;
    
    % Calcul de l'intégrale du terme source: int_K f
    intF = 0;
    aire=Aires(tri);
    [xyp,wp,lambda,np]=IntTri_Ham7(CoorNeuT);
    awp=aire.*wp';
    F = evalDonnees_Source(xyp(:,1),xyp(:,2));
    intF = sum(awp.*F ,1);

    int_f_values(tri) = intF;
    
    % Calcul de l'erreur d'équilibre
    erreur_equilibre(tri) = intF + int_gk_edge;
    
    % Détermination du type de triangle
    if ~isempty(find(RefEdg(TriEdg(tri,:)) > 0))
        triangle_type(tri) = 1; % Triangle sur bord
    else
        triangle_type(tri) = 0; % Triangle intérieur
    end
    
    % Vérification du critère DEBUG
    if abs(erreur_equilibre(tri)) > tolerance
        erreur_triangles = [erreur_triangles; tri];
        erreur_valeurs = [erreur_valeurs; erreur_equilibre(tri)];
    end
end

% Affichage des résultats
fprintf('Nombre total de triangles : %d\n', Nbtri);
fprintf('Nombre de triangles avec erreur > %.2e : %d\n', tolerance, length(erreur_triangles));
fprintf('Erreur maximale : %.6e\n', max(abs(erreur_equilibre)));
fprintf('Erreur moyenne : %.6e\n', mean(abs(erreur_equilibre)));

if ~isempty(erreur_triangles)
    fprintf('\nTriangles avec erreur :\n');
    for i = 1:length(erreur_triangles)
        tri = erreur_triangles(i);
        if triangle_type(tri) == 1
            type_str = 'Bord';
        else
            type_str = 'Intérieur';
        end
        fprintf('Triangle %d (%s): erreur = %.6e\n', tri, type_str, erreur_valeurs(i));
    end
end

%% VISUALISATION GRAPHIQUE
figure('Name', 'Analyse DEBUG - Triangles avec erreurs', 'Position', [100, 100, 1400, 800]);

% Graphique 1: Maillage avec triangles en erreur colorés
subplot(1,3,1);
triplot(NumTri, CoorNeu(:,1), CoorNeu(:,2), 'k-', 'LineWidth', 0.5);
hold on;

% Colorer les triangles en erreur
if ~isempty(erreur_triangles)
    for i = 1:length(erreur_triangles)
        tri = erreur_triangles(i);
        vertices = NumTri(tri, :);
        if triangle_type(tri) == 1
            fill(CoorNeu(vertices, 1), CoorNeu(vertices, 2), 'red', 'FaceAlpha', 0.7);
        else
            fill(CoorNeu(vertices, 1), CoorNeu(vertices, 2), 'blue', 'FaceAlpha', 0.7);
        end
        
        % Numéro du triangle au centre
        bary = mean(CoorNeu(vertices, :));
        text(bary(1), bary(2), num2str(tri), 'FontSize', 8, 'Color', 'white', 'FontWeight', 'bold');
    end
end

title('Triangles avec erreurs DEBUG');
xlabel('x'); ylabel('y');
legend('Maillage', 'Bord', 'Intérieur', 'Location', 'northeast');
axis equal; grid on;

% Graphique 2: Carte de chaleur des erreurs
subplot(1,3,2);
patch('Faces', NumTri, 'Vertices', CoorNeu, 'FaceVertexCData', log10(abs(erreur_equilibre) + 1e-16), ...
      'FaceColor', 'flat', 'EdgeColor', 'k', 'LineWidth', 0.3);
colorbar;
title('Log10(|Erreur équilibre|)');
xlabel('x'); ylabel('y');
axis equal;


% Graphique 4: int_gk vs int_f
subplot(1,3,3);
plot(int_f_values, -int_gk_values, 'bo', 'MarkerSize', 4);
hold on;
plot([min(int_f_values), max(int_f_values)], [min(int_f_values), max(int_f_values)], 'r--', 'LineWidth', 2);
xlabel('int_K f');
ylabel('-int_{∂K} gk');
title('Équilibre: int_K f + int_{∂K} gk = 0');
grid on;
axis equal;

% Marquer les triangles en erreur
if ~isempty(erreur_triangles)
    plot(int_f_values(erreur_triangles), -int_gk_values(erreur_triangles), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    legend('Tous triangles', 'Équilibre parfait', 'Triangles en erreur', 'Location', 'northeast');
end



endfunction 
