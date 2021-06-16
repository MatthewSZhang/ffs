function [IND, criteria] = ffs(X, Y, t, alg)
% Greedy Feature Selection by Orthogonal Least Squares
% Input:
% X: N * n matrix. N observations and n features. Feature matrix.
% Y: N * m matrix. N observations and m+1 classes. e.g. c-1 dummy encoded response.
% t: The number of features required to be selected.
% alg: choose h-correlation or theta-angle manually

% Output
IND = zeros(t, 1);
criteria = zeros(t, 1);

n = size(X, 2);
m = size(Y, 2);
INDr = 1:n; % index of rest features
N = size(Y, 1); % number of the observations

% Step 1
Xc = X - mean(X);
Yc = Y - mean(Y); 

% Step 2: Algorithm determination
if nargin < 4
    if N > m+n
        alg = 'theta';
    else
        alg = 'h';
    end
end

switch alg
    case 'h'
        Fx = Xc;
        Fy = Yc;
        Fn = N;
    case 'theta'
        [~, S, Rv] = svd([Xc, Yc], 'econ');
        F = S*Rv';
        Fx = F(:, 1:n);
        Fy = F(:, n+1:end);
        Fn = m + n;
end


% Step 3
V = f_orth(Fy);
Ws = zeros(Fn, t);
Wr = Fx./vecnorm(Fx);
for k = 1:t
    % Step 4
    g = Wr'*V;
    h = g.^2; % h-correlation or theta-angle
    RwV = sum(h, 2); % sum squared multiple correlation
    
    % Step 5
    [criteria(k), ind] = max(RwV);
    IND(k) = INDr(ind);
    INDr(:, ind) = [];

    % Step 3
    Ws(:, k) = Wr(:, ind);

    Wr(:, ind) = [];
    R = Ws(:, k)'*Wr;
    Wr = Wr - Ws(:, k).*R;
    Wr = Wr./vecnorm(Wr);
end

end


function [Q, R] = f_orth(X)
% Orthogonalisation of a Matrix by the Modified Gram-Schmidt
[N, n] = size(X);
Q = zeros(N, n);
R = zeros(n);
for k = 1:n
    Q(:, k) = X(:, k);
    for s = 1:k-1
        R(s, k) = Q(:, s)'*Q(:, k);
        Q(:, k) = Q(:, k) - R(s, k)*Q(:, s);
    end
    R(k, k) = norm(Q(:, k));
    Q(:, k) = Q(:, k)/R(k, k);
end
end


