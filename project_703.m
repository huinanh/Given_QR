%clear;
%clc;

% time
scale = [5 10 15 20 25 30 35 40 45 50];
naive = [0.016 0.046  0.117 0.167 0.319 0.491 1.520 1.631 2.116 3.780];
refined = [0.012 0.021 0.031 0.031 0.045 0.058 0.100 0.098 0.116 0.130];
global tol flag;
tol = 1e-8;
flag = 1;
% degrees of freedom
n = 40;

% for n = 35
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
%    A_ori = D;
%    error = 1;
%   count = 0;
%   flag = 1;
%   while flag
%        count = count+1;
%        flag = 0;
%        A = given_rotation_refined(A);
%    end
%    iteration(n-2) = count;
% end
% n = 3:20;
% plot(n,iteration)
n_time = [10 20 25 30 35 40]
time = [0.032 0.150 0.187 0.366 0.544 0.713]
time_shift = [0.015 0.046 0.053 0.093 0.152 0.209]
%iterate until A is close to upper triangular
for i = 1:1
%     disp("n = "+i)
%     A = given_rotation(A);
%     A = D;
%     A = given_rotation_refined(A);
    A = given_rotation_tri(A);
    A = given_rotation_shift(A);
end

%% naive given rotation
% conduct one iteration of QR decomposition and generate A for next
% iteration
function A = given_rotation(A)
    global flag tol;
    R = A;
   [m,n] = size(A);
    v = ones(1,m);
    Q_final = diag(v);
    for j = 1:n
        for i = m:-1:j+1
        % set elements small enough to zero and skip
            if abs(A(i,j)) < tol
                A(i,j) = 0;
                continue
            end
            flag = 1;
            [c,s] = givens(A(i-1,j),A(i,j));
            Q = Q_matrix(i-1,i,c,s,m);
            R = Q'*R;
            Q_final = Q_final * Q;
        end
    end
    A = R*Q_final;
end
%% refined given rotation

% conduct one iteration of QR decomposition and generate A for next
% iteration;
% refined: skip process when element is small enough;
%          only modify two rows under rotation instead of matrix
%          multiplication
function A = given_rotation_refined(A)
    global tol flag;
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
            A(i-1:i,j:n) =  [c s; -s c]'*A(i-1:i,j:n);
            %A(j:n,i-1:i) = A(j:n,i-1:i)*[c s;-s c];
            Q_final(1:n,i-1:i) = Q_final(1:n,i-1:i) * [c s;-s c];
        end
    end
    %Q_final'*Q_final;
    A = A*Q_final;
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
    index = 1;
    count = 1;
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
            %u = min(eig(A(n-1:n,n-1:n)));
            u = min(eig(A(1:2,1:2)));
            v = ones(1,n);
            I = diag(v);
            A = A - u*I;  % shift 
            [c,s] = givens(A(i-1,j),A(i,j));
            %Q = Q_matrix(i-1,i,c,s,m);
            A(i-1:i,j:n) =  [c s; -s c]'*A(i-1:i,j:n);
            A(j:n,i-1:i) = A(j:n,i-1:i)*[c s;-s c];
            Q_final(1:n,i-1:i) = Q_final(1:n,i-1:i) * [c s;-s c];
            %disp(count + " iteration");
            count = count + 1;
            A = A + u*I; % revert the shift
            end
        end
       
        if abs(A(2,1))<tol
            e(index) = A(1,1); % Ann converge to eigenvalue
            index = index + 1;
            eig(A);
            %A(n,n-1) = 0; A(n-1,n) = 0; A(n,n) = 0;
            %A = A(1:n-1,1:n-1); % reduce A down to n-1 X n-1
            A = A(2:n,2:n); % reduce A down to n-1 X n-1
        end
        [m,n] = size(A);
    end
    temp = eig(A);
    e(end - 1) = temp(1);
    e(end) = temp(2);
    e;
end



% generate c and s given the number rows that rotate, have different
% approaches
function [c,s] = givens(a,b)
    c = a/(sqrt(a^2+b^2));
    s = -1*b/(sqrt(a^2+b^2));
end

% generate rotation matrix Q
function Q = Q_matrix(i,j,c,s,m)
    v = ones(1,m);
    v(i) = c; v(j) = c;
    Q = diag(v);
    Q(i,j) = s;
    Q(j,i) = -1*s;
end
    
    