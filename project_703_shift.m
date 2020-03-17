%clear;
%clc;

global tol flag;
tol = 1e-6;



% degrees of freedom
n = 3;
%J_time = zeros(1,199);
%S_time = zeros(1,199);
naive_t = zeros(1,49);
refined_t = zeros(1,49);
for n = 2:50
flag = 1;

% mass and stiffness matrix
mi = ones(1,n);
M = diag(mi);
ki = ones(1,n);
K = zeros(n,n);

K(1,1) = ki(1)+ki(2);
K(2,1) = -ki(2);
for i=2:n-1
   K(i-1,i) = -ki(i);
   K(i,i) = ki(i)+ki(i+1);
   K(i+1,i) = -ki(i+1);
end
K(n-1,n) = -ki(n);
K(n,n) = ki(n);

D = K^-1*M;
A = D;

% iterate till converge
t = cputime;
while flag
    A = given_rotation(A);
end
e = cputime - t;
naive_t(n-1) = e;

t = cputime;
flag = 1;
A = D;
while flag
    A = given_rotation_refined(A);
end
e = cputime - t;
refined_t(n-1) = e;
end
% iterate until A is close to upper triangular
% t = cputime;
% A = given_rotation_tri(A);
% A = given_rotation_shift(A);
% e = cputime-t;
% S_time(n-1) = e;
% end

n = 2:200;


%% naive given rotation
% conduct one iteration of QR decomposition and generate A for next
% iteration
function A = given_rotation(A)
    global tol flag;
    flag = 0;
    R = A;
   [m,n] = size(A);
    v = ones(1,m);
    Q_final = diag(v);
    for j = 1:n
        for i = m:-1:j+1
            if abs(A(i,j)) > tol
                flag = 1;
            end
            
            [c,s] = givens(A(i-1,j),A(i,j));
            Q = Q_matrix(i-1,i,c,s,m);
            A = Q'*A*Q;
            Q_final = Q_final * Q;
        end
    end
    %A = R*Q_final;
end
%% refined given rotation

% conduct one iteration of QR decomposition and generate A for next
% iteration;
% refined: skip process when element is small enough;
%          only modify two rows under rotation instead of matrix
%          multiplication
function A = given_rotation_refined(A)
    global tol flag;
    flag = 0;
   [m,n] = size(A);
    v = ones(1,m);
    Q_final = diag(v);
    for j = 1:n-1
        for i = m:-1:j+1
            % set elements small enough to zero and skip
            if abs(A(i,j)) < tol
                A(i,j) = 0;
                continue
            end
            flag = 1;
            [c,s] = givens(A(i-1,j),A(i,j));
            %Q = Q_matrix(i-1,i,c,s,m);
            %A(i-1:i,j:n) =  [c s; -s c]'*A(i-1:i,j:n)
            %A(j:n,i-1:i) = A(j:n,i-1:i)*[c s;-s c]
            A(i-1:i,1:n) =  [c s; -s c]'*A(i-1:i,1:n);
            A(1:n,i-1:i) = A(1:n,i-1:i)*[c s;-s c];
            %Q_final(1:n,i-1:i) = Q_final(1:n,i-1:i) * [c s;-s c];
        end
    end
end

%% given rotation that reduces matrix down to tridiagonal form
% conduct one iteration of QR decomposition and generate A for next
% iteration;
% refined: skip process when element is small enough;
%          only modify two rows under rotation instead of matrix
%          multiplication
% reduce symmetrxi matrix to tridiagonal form
function A = given_rotation_tri(A)
    global tol;
   [m,n] = size(A);
    v = ones(1,m);
    Q_final = diag(v);
    for j = 1:n-2
        for i = m:-1:j+2
        % set elements small enough to zero and skip
        if abs(A(i,j)) < tol
            A(i,j) = 0;
            continue
        end
        [c,s] = givens(A(i-1,j),A(i,j));
        %Q = Q_matrix(i-1,i,c,s,m);
        A(i-1:i,j:n) =  [c s; -s c]'*A(i-1:i,j:n);
        A(j:n,i-1:i) = A(j:n,i-1:i)*[c s;-s c];
        Q_final(1:n,i-1:i) = Q_final(1:n,i-1:i) * [c s;-s c];
        end
    end
end



%% given rotation with shift
% conduct one iteration of QR decomposition and generate A for next
% iteration;
% refined: skip process when element is small enough;
%          only modify two rows under rotation instead of matrix
%          multiplication
% reduce symmetrxi matrix to tridiagonal form
function A = given_rotation_shift(A)
    global tol;
   [m,n] = size(A);
    e = ones(1,n); % store all eigenvalues
    v = ones(1,m);
    Q_final = diag(v); % initial Q equal to identical matrix
    % have a guess (shift as eigenvalue)
    while(n>2)
        
        % one iteration of QR
        for j = 1:n
            for i = m:-1:j+1
            % set elements small enough to zero and skip
            if abs(A(i,j)) < tol
                A(i,j) = 0;
                continue
            end
            u = min(eig(A(n-1:n,n-1:n)));
            v = ones(1,n);
            I = diag(v);
            A = A - u*I;  % shift 
            [c,s] = givens(A(i-1,j),A(i,j));
            %Q = Q_matrix(i-1,i,c,s,m);
            A(i-1:i,1:n) =  [c s; -s c]'*A(i-1:i,1:n);
            A(1:n,i-1:i) = A(1:n,i-1:i)*[c s;-s c];
            Q_final(1:n,i-1:i) = Q_final(1:n,i-1:i) * [c s;-s c];
            A = A + u*I; % revert the shift
            end
        end
       
        if abs(A(n,n-1))<tol
            e(n) = A(n,n); % Ann converge to eigenvalue
            %eig(A)
            %A(n,n-1) = 0; A(n-1,n) = 0; A(n,n) = 0;
            A = A(1:n-1,1:n-1); % reduce A down to n-1 X n-1
        end
        [m,n] = size(A);
    end
    temp = eig(A);
    e(1) = temp(1);
    e(2) = temp(2);
    e;
end



% generate c and s given the number rows that rotate, have different
% approaches
function [c,s] = givens(a,b)
    c = -1*a/(sqrt(a^2+b^2));
    s = b/(sqrt(a^2+b^2));
end

% generate rotation matrix Q
function Q = Q_matrix(i,j,c,s,m)
    v = ones(1,m);
    v(i) = c; v(j) = c;
    Q = diag(v);
    Q(i,j) = s;
    Q(j,i) = -1*s;
end
    
    