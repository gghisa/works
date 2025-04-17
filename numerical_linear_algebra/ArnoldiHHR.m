function [V,H,L] = ArnoldiHHR(A,v,k)

[n, ~] = size(A);
L = [];
V = [];
H = [];
u = v;

% now continue looking for l2

for j=1:k+1

    e_j = [zeros(j - 1,1); ones(1); zeros(n - j,1)];

    % construct u_hat and l

    u_hat = u(j:n);

    sign = 1;
    if u_hat(1)<0
        sign = -1;
    end

    h = [ u(1:j - 1); - sign * norm(u_hat); zeros(n - j, 1)];
    
    if j>1
        H = [H h];
    end

    l  = u_hat + sign * norm(u_hat) * [ 1; zeros(n - j,1) ];

    l = [zeros(j - 1,1); l];

    L = [L l]; % from here on we use l from L

    % find v_j

    y = e_j - 2 * L(:,j) * ( L(j,j) / ( L(:,j)' * L(:,j) ) );

    for i=1:j-1
        y = y - 2 * L(:,j - i) * ( L(:,j - i)' * y) / ( L(:,j - i)' * L(:,j - i));
    end
    % at the end of the above loop, y is one of our basis vectors

    V = [V y]; % from here on we use y from V
    
    % we now produce the next u
    u = A * y;

    for i=1:j
        u = u - 2 * L(:,i) * (L(:,i)' * u) / (L(:,i)' * L(:,i));
    end
end