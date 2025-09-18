function coeffs = robust_null_combined(M)
    [U, S, V] = svd(M);
    s = diag(S);
    n = size(M, 2);
    
    % Critères multiples pour identifier le noyau
    tol_abs = 1e-10;
    tol_rel = max(size(M)) * eps(max(s));
    tol = max(tol_abs, tol_rel);
    
    % Trouver le nombre de valeurs singulières significatives
    rank_M = sum(s > tol);
    
    if rank_M < n
        % Il y a un noyau non-trivial
        null_dim = n - rank_M;
        if null_dim == 1
            coeffs = V(:, end);
        else
            % Plusieurs dimensions dans le noyau, prendre une combinaison
            % ou la première dimension
            coeffs = V(:, rank_M + 1);
        end
    else
        % Rang numériquement plein, prendre le vecteur correspondant 
        % à la plus petite valeur singulière
        coeffs = V(:, end);
    end
    
    % S'assurer que la solution n'est pas trop petite
    if norm(coeffs) < 1e-14
        % Solution de secours : vecteur canonique
        coeffs = zeros(n, 1);
        coeffs(1) = 1;
    else
        coeffs = coeffs / norm(coeffs);
    end
    
    % Vérification finale : s'assurer que M*coeffs est petit
    residual = norm(M * coeffs);
    if residual > 1e-6
        warning('Solution trouvée avec résidu élevé: %e', residual);
    end
end

