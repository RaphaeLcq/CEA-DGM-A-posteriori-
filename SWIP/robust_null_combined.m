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
% Author : Raphael Lecoq, CEA
%
%
% SYNOPSIS
% Robust method to find the kernel of a given matrix.
%
% GLOBAL : None
% INPUT :
%   - M        : matrix
% OUTPUT:
%   - coeffs    : coefficients in kernel(M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function coeffs = robust_null_combined(M)
    [U, S, V] = svd(M);
    s = diag(S);
    n = size(M, 2);

    % Different criteria to find the kernel
    tol_abs = 1e-10;
    tol_rel = max(size(M)) * eps(max(s));
    tol = max(tol_abs, tol_rel);

    % Find the significative roots
    rank_M = sum(s > tol);

    if rank_M < n
        % Theres a non trivial kernel
        null_dim = n - rank_M;
        if null_dim == 1
            coeffs = V(:, end);
        else
            % The kernel is more than one dimension thus take the last rank
            % or the first rank
            coeffs = V(:, rank_M + 1);
        end
    else
        % If the matrix is fullrank then take the smallest eingenvalue
        coeffs = V(:, end);
    end

    % Sanity check over the norm coeffs
    if norm(coeffs) < 1e-14
        % If the norm coeffs are 0, canonical vector
        coeffs = zeros(n, 1);
        coeffs(1) = 1;
    else
        coeffs = coeffs / norm(coeffs);
    end

    % Sanity check ; norm(Mx) < tol
    residual = norm(M * coeffs);
    if residual > 1e-6
        warning('Solution trouvée avec résidu élevé: %e', residual);
    end
end

