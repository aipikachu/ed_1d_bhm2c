% MATLAB script illustrating the use of MEXPV.

disp('Loading the matrix ...');
A = loadcrs('../data/c1024.crs');

disp('Computing w = exp(A)e_1 ...');
[n,n] = size(A);
v = eye(n,1);
[w,err] = mexpv(1,A,v);

disp('w(1:10) =');
disp(w(1:10));

disp('err =');
disp(err)


% generate a transition rate matrix
n = 20000;
A = rand(n);
for j = 1:n
    sumj = 0;
    for i = 1:n
        if rand < 0.5, A(i,j) = 0; end;
        sumj = sumj + A(i,j);
    end;
    A(j,j) = A(j,j)-sumj;
end;
v = eye(n,1);
A = sparse(A); % invaluable for a large and sparse matrix.

tic
[w,err] = expv(1,A,v);
toc

disp('w(1:10) ='); disp(w(1:10));
disp('err =');     disp(err);

% tic
% w_matlab = expm(full(A))*v;
% toc
% 
% disp('w_matlab(1:10) ='); disp(w_matlab(1:10));
% gap = norm(w-w_matlab)/norm(w_matlab);
% disp('||w-w_matlab|| / ||w_matlab|| ='); disp(gap);

% spy(A)