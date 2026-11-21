function Ak = rank_k_approx(U, s, V, k)
%RANK_K_APPROX  Best rank-k approximation of a matrix from its SVD.
%   Ak = RANK_K_APPROX(U, s, V, k) returns
%       Ak = sum_{i=1..k} s(i) * U(:,i) * V(:,i)'
%   where [U, S, V] = svd(A, 'econ') and s = diag(S).
%   By the Eckart-Young theorem, no other rank-k matrix is closer to A
%   in the Frobenius norm.

    Ak = U(:, 1:k) * diag(s(1:k)) * V(:, 1:k)';
end
