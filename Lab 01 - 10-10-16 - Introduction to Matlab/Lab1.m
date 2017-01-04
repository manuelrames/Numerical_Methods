clear
clc

%% Build matrix A using a for loop

n = 10;
h = 1/n;
dim = n-1;
A = zeros(dim,dim);

for i = 1:dim
    A(i,i) = 2;
end

for i = 1:dim-1
    A(i,i+1) = -1;
    A(i+1,i) = -1;
end

A = h^(-2) * A;

%% Build matrix A using command diag
clear
clc

n = 10;
h = 1/n;
A = 2*diag(ones(n-1,1))-diag(ones(n-2,1),-1)-diag(ones(n-2,1),1);
A = h^(-2) * A;
whos A % see memory consumption of A
A_sparse = sparse(A);
whos A_sparse

%% Build matrix A in sparse format
clear
clc

n = 10;
h = 1/n;
e = ones(n-1,1);
A_sparse = h^(-2) * spdiags([2*e,-e,-e],[0,-1,1],n-1,n-1);
%A = full(A_sparse) % to go back to the full format
spy(A_sparse) % to see the pattern of A
nnz(A_sparse) % to see how many non-zero elements there are in A