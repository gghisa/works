function [V,T,Ys,MU] = Lanczos_EIGS(A,v)

[n,~] = size(A);

V = v/norm(v); % initialise and normalise the basis matrix

w = A*V(:,1); % apply A to the first basis vector

% project A*v on v
T(1,1) = V(:,1)'*w; % store the projection coefficient
w = w - V(:,1)*T(1,1); % orthogonalise A*v wrt v
T(2,1) = norm(w); % store the norm of orthogonalised A*V
w = w/T(2,1); % normalise it
V = [V w]; % append second orthonormal basis vector

% eigenproblem
[Ys, MU] = eig(T(1,1)); % get eigenvectors and eigenvalues for T_k,k
P = [1;T(2,1)]; % initialise print matrix
r_1_norm = abs(P(2,1)) * Ys(end,1);
P = [P, [MU(1,1); r_1_norm]];
P(:,2:end) = sortrows(P(:,2:end)',1)';
P
% plot
scatter(1,P(1,2),'x','blue','LineWidth',3)
title('Convergence of Ritz values, made by Giovanni and Svenja')
xlabel('Iteration #')
hold on

for j=2:n
    T(j-1,j) = T(j,j-1); % could be V(:,j-1)'*w but we want exact T symmetry
    w = A*V(:,j);
    T(j,j) = V(:,j)'*w;
    % apply MGS orthogonalisation
    for i=j-1:j
        w = w - V(:,i) * (V(:,i)' * w);
    end
    T(j+1,j) = norm(w);
    w = w/T(j+1,j);
    V = [V w];
    % eigenproblem
    [Ys, MU] = eig(T(1:j,1:j)); % get eigenvectors and eigenvalues for T_k,k
    P = [j;T(j+1,j)]; % initialise print matrix
    for q=1:j
        r_q_norm = abs(P(2,1) * Ys(end,q));
        P = [P, [MU(q,q); r_q_norm]];
    end
    P(:,2:end) = sortrows(P(:,2:end)',1)';
    P
    % plot
    scatter(j,P(1,2:end),'x','red','LineWidth',3)    
end
end