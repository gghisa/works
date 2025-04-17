n = 10;

%A = gallery('grcar',n);
%v = ones(n,1);

A = diag(0:10);
v = (2.^(0:10))';

V = []; H = []; L = [];

running = true;

while running
    
    if isempty(V)
        disp('Empty arnoldi basis, would you like to extend the arnoldi basis? [Y/N]');
        x = input('', 's');
        if strcmp(x,'Y')
            [V,H] = ExtendArnoldi(A,v/norm(v),H);
            disp('Current Ritz data:');
            L = ListRitzData(H);
            disp(L);
        end
    else
        disp('Would you like to EXTEND, FILTER, or EXIT?')
        x = input('', 's');

        if strcmp(x,'EXTEND')
            [V,H] = ExtendArnoldi(A,V,H);
            disp('Current Ritz data:');
            L = ListRitzData(H);
            disp(L);

        elseif strcmp(x,'FILTER')
            disp('Enter eigenvalue index for that which you wish to remove.')
            disp("- Indeces are matched to Ritz values' display")
            disp('- Enter -1 if you wish to filter no more eigenvalues.')
            y = 0;
            idxs = [];
            while y >= 0
                y = input('Enter (integer) index: ');
                if y > 0
                    idxs = [idxs, y];
                end
            end
            for i=1:length(idxs)
                disp("Removing eigenvalue " + idxs(i) + ": " + L(1,1+idxs(i)))
                [V,H] = FilterAway(L(1,1+idxs(i)), V, H);
                disp('Current Ritz data:');
                L = ListRitzData(H);
                disp(L);
            end

        elseif strcmp(x,'EXIT')
            running = false;
        end
    end
end