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
% SYNOPSIS Affiche la solution calculée
% GLOBALS:
%        - CoorNeu(Nbpt,2) : coordonnees (x, y) des sommets (noeuds P1)
%        - Nbtri  : nombre de triangle
%        - CoorBary(Nbtri,2)   : Coordonnees des barycentres de triangles
%        - NumTri(Nbtri,3) : liste de triangles  (3 numeros de sommets) 
%        - Nbtri  : nombre de triangle
%        - Phih(idof, 3*ordre) : solution discrete dans la base dG

function val = evalSolutionDiscrete(x, y)


    global CoorNeu
    global Nbtri CoorBary NumTri
    global Phih

    % S'assurer que x et y ont la même taille
    if isscalar(x) && ~isscalar(y)
        x = x * ones(size(y));
    elseif ~isscalar(x) && isscalar(y)
        y = y * ones(size(x));
    end

    % Convertir en vecteurs colonnes pour traitement
    original_size = size(x);
    x_vec = x(:);
    y_vec = y(:);
    n_points = length(x_vec);

    % Initialiser le résultat
    val_vec = zeros(n_points, 1);

    % Version optimisée : pré-calculer les données des triangles
    if ~exist('triangle_data_precomputed', 'var')
        triangle_data = precompute_triangle_data();
    else
        global triangle_data_precomputed;
        triangle_data = triangle_data_precomputed;
    end

    % Traitement vectorisé par chunks pour éviter les problèmes de mémoire
    chunk_size = min(1000, n_points);  % Traiter par chunks de 1000 points max

    for chunk_start = 1:chunk_size:n_points
        chunk_end = min(chunk_start + chunk_size - 1, n_points);
        chunk_indices = chunk_start:chunk_end;

        x_chunk = x_vec(chunk_indices);
        y_chunk = y_vec(chunk_indices);

        % Trouver les triangles pour ce chunk
        triangle_indices = find_triangles_vectorized(x_chunk, y_chunk, triangle_data);

        % Évaluer la solution pour chaque point du chunk
        for i = 1:length(chunk_indices)
            if triangle_indices(i) > 0
                val_vec(chunk_indices(i)) = eval_point_in_triangle(...
                    x_chunk(i), y_chunk(i), triangle_indices(i));
            else
                val_vec(chunk_indices(i)) = 0;  % ou NaN selon votre préférence
            end
        end
    end

    % Remodeler selon la forme originale
    val = reshape(val_vec, original_size);
end

function triangle_data = precompute_triangle_data()
    % Pré-calcule les données géométriques des triangles
    global CoorNeu NumTri Nbtri

    triangle_data = struct();
    triangle_data.x1 = CoorNeu(NumTri(:,1), 1);
    triangle_data.y1 = CoorNeu(NumTri(:,1), 2);
    triangle_data.x2 = CoorNeu(NumTri(:,2), 1);
    triangle_data.y2 = CoorNeu(NumTri(:,2), 2);
    triangle_data.x3 = CoorNeu(NumTri(:,3), 1);
    triangle_data.y3 = CoorNeu(NumTri(:,3), 2);

    % Pré-calculer les dénominateurs
    triangle_data.denom = (triangle_data.y2 - triangle_data.y3) .* ...
                         (triangle_data.x1 - triangle_data.x3) + ...
                         (triangle_data.x3 - triangle_data.x2) .* ...
                         (triangle_data.y1 - triangle_data.y3);

    % Marquer les triangles valides (non dégénérés)
    triangle_data.valid = abs(triangle_data.denom) >= 1e-12;
end

function triangle_indices = find_triangles_vectorized(x_vec, y_vec, triangle_data)
    % Trouve les triangles contenant les points de manière vectorisée
    global Nbtri

    n_points = length(x_vec);
    triangle_indices = zeros(n_points, 1);
    tolerance = 1e-10;

    % Pour chaque point, tester tous les triangles
    for i = 1:n_points
        x = x_vec(i);
        y = y_vec(i);

        % Calcul vectorisé des coordonnées barycentriques pour tous les triangles
        valid_triangles = find(triangle_data.valid);

        if isempty(valid_triangles)
            continue;
        end

        % Extraire les données pour les triangles valides
        x1 = triangle_data.x1(valid_triangles);
        y1 = triangle_data.y1(valid_triangles);
        x2 = triangle_data.x2(valid_triangles);
        y2 = triangle_data.y2(valid_triangles);
        x3 = triangle_data.x3(valid_triangles);
        y3 = triangle_data.y3(valid_triangles);
        denom = triangle_data.denom(valid_triangles);

        % Calcul vectorisé des coordonnées barycentriques
        lambda1 = ((y2 - y3) .* (x - x3) + (x3 - x2) .* (y - y3)) ./ denom;
        lambda2 = ((y3 - y1) .* (x - x3) + (x1 - x3) .* (y - y3)) ./ denom;
        lambda3 = 1 - lambda1 - lambda2;

        % Vérifier quels triangles contiennent le point
        inside = (lambda1 >= -tolerance) & (lambda2 >= -tolerance) & (lambda3 >= -tolerance);

        % Prendre le premier triangle trouvé
        triangle_found = find(inside, 1);
        if ~isempty(triangle_found)
            triangle_indices(i) = valid_triangles(triangle_found);
        end
    end
end

function val = eval_point_in_triangle(x, y, triangle_idx)
    % Évalue la solution en un point d'un triangle donné
    global Phih idof

    debT = idof(triangle_idx, 1);
    finT = idof(triangle_idx, 2);
    Uh_T = Phih(debT:finT);

    % Évaluer les fonctions de base
    base_values = eval_fct_base([x, y], triangle_idx);
    val = base_values' * Uh_T;
end
