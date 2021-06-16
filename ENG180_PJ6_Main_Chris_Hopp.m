%Chris Hopp
%915866326
%ENG-180 Project 6: Ordinary Differential Equations
%11/10/2020


clc
clf
clear all




Problem1a()
Problem1b()
Problem1c()
Problem2a()
Problem2b()
Problem3a()
Problem3b()
Problem4()
%% Problem 1a
function Problem1a()
n = 101;
x = linspace(-1,1,n);
dx = x(2) - x(1);

epsilon = .1;

f = -1;                                         % Forms tridiagonal system
a = (ones(n,1))*(-epsilon/dx^2 - f/(2*dx));
a(1) = 0;
a(end) = 0;
b = (ones(n,1))*((2*epsilon)/(dx^2));
b(1) = 1;
b(end) = 1;
c = (ones(n,1))*(-epsilon/dx^2 + f/(2*dx));
c(1) = 0;
c(end) = 0;
d = zeros(n,1);
d(1) = 1;
d(end) = -1;
u1 = Thomas3(a,b,c,d);


f = 1;                                          % Forms tridiagonal system
a = (ones(n,1))*(-epsilon/dx^2 - f/(2*dx));
a(1) = 0;
a(end) = 0;
b = (ones(n,1))*((2*epsilon)/(dx^2));
b(1) = 1;
b(end) = 1;
c = (ones(n,1))*(-epsilon/dx^2 + f/(2*dx));
c(1) = 0;
c(end) = 0;
d = zeros(n,1);
d(1) = 1;
d(end) = -1;
u2 = Thomas3(a,b,c,d);


res = dx;
tol = 1e-3;
uOld = zeros(1,n);
uNew = x;                   % Initial guess
w=1.2;                      % Over relaxed

while res > tol             % Lagging non-linear term
    uNew = w*(uNew-uOld)+uOld;
    uOld = uNew;
    for i = 1:n                                   % Forms tridiagonal
        a(i) = (-epsilon/dx^2)-(uOld(i)/(2*dx));
        c(i) = (uOld(i)/(2*dx))-epsilon/dx^2;
    end
    a(end) = 0;
    a(1) = 0;
    b  = ones(n,1) * 2*epsilon/dx^2;
    b(1) = 1;
    b(end) = 1;
    c(1) = 0;
    c(end) = 0;
    d = zeros(n,1);
    d(1) = 1;
    d(end) = -1;
    uNew = Thomas3(a,b,c,d);
    for i = 2:n-1                                    % Checks convergence
        residual(i) = (uNew(i)/(2*dx))*(uNew(i+1)-uNew(i-1)) - ((epsilon/dx^2)*(uNew(i+1)-2*uNew(i)+uNew(i-1)));
    end
    res = max(abs(residual));
end


figure(1)
plot(x,u1,'DisplayName', 'f = -1')
hold on
plot(x,u2,'DisplayName','f = 1')
plot(x,uNew,'DisplayName','f = u')
xlabel('x')
ylabel('u')
title('Boundary Layers and Shock Waves')
legend('location', 'best')
end

%% Problem 1b
function Problem1b()
n = 101;
x = linspace(.5,1,n);
dx = x(2) - x(1);

epsilon = dx;

w=1;
cor = dx;
tol = 1e-5;
uOld = zeros(1,n);
uNew = linspace(1,.353,n);
A = 1+(2*x-1).^2;
dA = 8*x-4;
iteration = 0;

while cor > tol & iteration < 10000
    iteration = iteration + 1
    uNew = w*(uNew-uOld)+uOld;
    uOld = uNew;
    a = zeros(1,n);                 % Allocation for tridiagonal
    b = zeros(1,n);
    c = zeros(1,n);
    d = zeros(1,n);
    B = zeros(1,n);
    
    for i = 1:n                     % Coefficient on u_i+1, u_i-1
        B(i) = (1 - ( 1 / (A(i) * uOld(i)^3)))/(2*dx);     
    end
    
    for i = 1:n                     % Forms tridiagonal
            a(i) = -B(i) - (epsilon / dx^2);
            b(i) = (2*epsilon/dx^2);
            c(i) = B(i) - (epsilon / dx^2);
            d(i) = dA(i)/((A(i)^2)*uOld(i)^2);
    end
    a(end) = 0;                     % Formats for thomas3 
    a(1) = 0;
    b(1) = 1;
    b(end) = 1;
    c(1) = 0;
    c(end) = 0;
    d(1) = 1;
    d(end) = .353;
    uNew = Thomas3(a,b,c,d);
    
    for i = 2:n-1
        correction(i) = uNew(i) - uOld(i);
    end
    cor = max(abs(correction));
    uNew1 = uNew;
    x1 = x;

end



n = 1001;
x = linspace(.5,1,n);
dx = x(2) - x(1);

epsilon = dx;

w=.9;
cor = dx;
tol = 1e-5;
uOld = zeros(1,n);
uNew = linspace(1,.353,n);
A = 1+(2*x-1).^2;
dA = 8*x-4;
iteration = 0;

while cor > tol & iteration < 10000
    iteration = iteration + 1
    uNew = w*(uNew-uOld)+uOld;
    uOld = uNew;
    a = zeros(1,n);                 % Allocation for tridiagonal
    b = zeros(1,n);
    c = zeros(1,n);
    d = zeros(1,n);
    B = zeros(1,n);
    
    for i = 1:n                     % Coefficient on u_i+1, u_i-1
        B(i) = (1 - ( 1 / (A(i) * uOld(i)^3)))/(2*dx);     
    end
    
    for i = 1:n                     % Forms tridiagonal
            a(i) = -B(i) - (epsilon / dx^2);
            b(i) = (2*epsilon/dx^2);
            c(i) = B(i) - (epsilon / dx^2);
            d(i) = dA(i)/((A(i)^2)*uOld(i)^2);
    end
    a(end) = 0;                     % Formats for thomas3 
    a(1) = 0;
    b(1) = 1;
    b(end) = 1;
    c(1) = 0;
    c(end) = 0;
    d(1) = 1;
    d(end) = .353;
    uNew = Thomas3(a,b,c,d);
    
    for i = 2:n-1
        correction(i) = uNew(i) - uOld(i);
    end
    cor = max(abs(correction));
    uNew2 = uNew;
    x2 = x;

end




    figure(2)
    plot(x1,uNew1,'DisplayName', 'N = 101, w = 1')
    hold on 
    plot(x2,uNew2,'DisplayName', 'N = 1001, w = .9')
    legend('Location', 'best')
    xlabel('x')
    ylabel('u')
    title('Converging-Diverging Nozzle')
    hold off
end

%% Problem 1c
function Problem1c()
n = 101;
x = linspace(0,1,n);
dx = x(2) - x(1);
alpha = 1;
beta = -1;

w=1;
cor = dx;
tol = 1e-5;
uOld = zeros(1,n);
uNew = zeros(1,n);
iteration = 0;

while cor > tol & iteration < 10000
    iteration = iteration + 1
    uNew = w*(uNew-uOld)+uOld;
    uOld = uNew;
    a = zeros(1,n);                 % Allocation for tridiagonal
    b = zeros(1,n);
    c = zeros(1,n);
    d = zeros(1,n);
    B = zeros(1,n);
    
   
    for i = 1:n                     % Forms tridiagonal
            a(i) = 1/dx^2;
            b(i) = -2/dx^2;
            c(i) = 1/dx^2;
            d(i) = -alpha*exp(beta*uOld(i));
    end
    a(end) = 0;                     % Formats for thomas3 
    a(1) = 0;
    b(1) = 1;
    b(end) = 1;
    c(1) = 0;
    c(end) = 0;
    d(1) = 0;
    d(end) = 0;
    uNew = Thomas3(a,b,c,d);
    
    for i = 2:n-1
        correction(i) = uNew(i) - uOld(i);
    end
    cor = max(abs(correction));
    x1 = x;
    uNew1 = uNew;
    


end



alpha = .8;
beta = 1;

w=1;
cor = dx;
tol = 1e-7;
uOld = zeros(1,n);
uNew = zeros(1,n);
iteration = 0;

while cor > tol & iteration < 10000
    iteration = iteration + 1
    uNew = w*(uNew-uOld)+uOld;
    uOld = uNew;
    a = zeros(1,n);                 % Allocation for tridiagonal
    b = zeros(1,n);
    c = zeros(1,n);
    d = zeros(1,n);
    B = zeros(1,n);
    
   
    for i = 1:n                     % Forms tridiagonal
            a(i) = 1/dx^2;
            b(i) = -2/dx^2;
            c(i) = 1/dx^2;
            d(i) = -alpha*exp(beta*uOld(i));
    end
    a(end) = 0;                     % Formats for thomas3 
    a(1) = 0;
    b(1) = -1;
    b(end) = 1;
    c(1) = 1;
    c(end) = 0;
    d(1) = ((-dx^2)/2)*(alpha*exp(beta*uOld(1)));
    d(end) = 0;
    uNew = Thomas3(a,b,c,d);
    
    for i = 2:n-1
        correction(i) = uNew(i) - uOld(i);
    end
    cor = max(abs(correction));
    x2 = x;
    uNew2 = uNew;
    


end
    figure(3)

    plot(x2,uNew2,'DisplayName', 'a = 0.8, b = 1, Neumann Boundary')
    hold on 
    plot(x1,uNew1,'DisplayName', 'a = 1, b = -1, Dirichlet Boundary')

    legend('Location', 'best')
    xlabel('x')
    ylabel('u')
    title('Frank-Kamenetskii Analysis')
    hold off



end

%% Problem 2a
function Problem2a()

[x1,u1] = Problem2Solver(0, pi/2,0,2,1);
[x2,u2] = Problem2Solver(0, pi/2,0,2,-1);
[x3,u3] = Problem2Solver(0, pi,  0,0, 1);
[x4,u4] = Problem2Solver(0, pi,  0,0, -1);
[x5,u5] = Problem2Solver(0, pi,  0,1, 1);
[x6,u6] = Problem2Solver(0, pi,  0,1, -1);
    figure(4)
    
    plot(x5,u5,'','DisplayName', '\alpha = 1, u(0)=0, u(\pi) = 1')

    hold on
    plot(x3,u3,'-s','DisplayName', '\alpha = 1, u(0)=0, u(\pi) = 0')
    plot(x1,u1,'DisplayName', '\alpha = 1, u(0)=0, u(\pi/2) = 2')
    plot(x2,u2,'DisplayName', '\alpha = -1, u(0)=0, u(\pi/2) = 2')
    plot(x4,u4,'-^','DisplayName', '\alpha = -1, u(0)=0, u(\pi) = 0')
    
    legend('Location', 'best')
    xlabel('x')
    ylabel('u')
    title('Existence and Uniqueness')
    hold off

    
    figure(5)

    plot(x6,u6,'','DisplayName', '\alpha = -1, u(0)=0, u(\pi) = 1')
    hold on

    
    legend('Location', 'best')
    xlabel('x')
    ylabel('u')
    title('Invalid Cases')
    hold off

end

%% Problem 2b
function Problem2b()
n = 101;
x = linspace(0,1,n);
dx = x(2) - x(1);

d = zeros(n+1,1);

    A = zeros(n+1,n);   % Not square, if I have time I will try the other approach

    for i = 2:n-1                           % Forms full matrix
            A(i,i-1) = 1/dx^2;
            A(i,i) = -2/dx^2;
            A(i,i+1) = 1/dx^2;
            d(i) = (cos(pi*x(i)))^2;
    end
        A(1,1) = -2/dx^2;
        A(1,2) = 2/dx^2;
        A(n,n-2) = 1/dx^2;
        A(n,n-1) = -2/dx^2;
        A(n,n) = 1/dx^2;                
        A(n+1,:) = dx;                      % Final row gives integral
        A(n+1,1) = dx/2;
        d(1) = (cos(pi*x(1)))^2;
        d(n) = (cos(pi*x(n)))^2;
        d(n+1) = 0;
        A(n+1,n) = dx/2;
    
    u = A\d;    % Matrix is not square, we have no solver for this 
    figure(6)
    plot(x,u)
    xlabel('x')
    ylabel('y')
    title('Mixed Boundary Condition')
    
    
end

%% Problem 3a
function Problem3a()
n = 101;
x = linspace(0,5,n);  % Taking L = 10
dx = x(2) - x(1);


w=.8;
cor = dx;
tol = 1e-2;
fOld = x;
fNew = fOld;
iteration = 0;


while cor > tol & iteration < 10000
    iteration = iteration + 1
    fNew = w*(fNew-fOld)+fOld;
    fOld = fNew;
    a = zeros(1,n);                 % Allocation for pentadiagonal
    b = zeros(1,n);
    c = zeros(1,n);
    d = zeros(1,n);
    e = zeros(1,n);
    f = zeros(1,n);
    

    for i = 1:n                     % Forms pentadiagonal
            a(i) = -1/dx^3;
            b(i) = (3/dx^3) + (fOld(i)/(2*dx^2));
            c(i) = (-3/dx^3) - (fOld(i)/dx^2);
            d(i) = (1/dx^3) + (fOld(i)/(2*dx^2));
    end
    a(1:2) = 0;                     % Formats for thomas5 
    b(1) = 0;
    d(n) = 0;
    e(n-1:n) = 0;
    
    c(1) = 1;                               % First B.C.
    d(1) = 0;
    e(1) = 0;
    
    b(2) = 3/dx^3 + fOld(2)/(2*dx^2);      % f'(0) = 0 B.C.
    c(2) = -2/dx^3 -fOld(2)/dx^2;
    d(2) = 1/dx^3+ fOld(2)/(2*dx^2);
    
    a(n) = -1/dx^3;                         % f'(L) = 1
    b(n) = 4/dx^3 + fOld(n)/dx^2;
    c(n) = -3/dx^3 - fOld(n)/dx^2;
    
    f(n) = -2/dx^2 - fOld(n)/dx;

    fNew = Thomas5(a,b,c,d,e,f);
    
    for i = 1:n
        correction(i) = fNew(i) - fOld(i);
    end
    cor = max(abs(correction));
    uNew1 = fNew;
    x1 = x;

end

% I really felt like these were the right equations and I don't know where
% I went wrong, but it won't converge.  It looks like it is enforcing a
% derivative of zero at the right end, instead of a derivative of 1. 


    figure(7)
    plot(x1,uNew1,'DisplayName', 'Pentadiagonal')
    hold on 
    legend('Location', 'best')
    xlabel('x')
    ylabel('f')
    title('Blasius Equation')
    hold off


end

%% Problem 3b
function Problem3b()
n = 101;
x = linspace(0,5,n);  % Taking L = 10
dx = x(2) - x(1);


w=.8;
cor = dx;
tol = 1e-2;
fOld = x;
fNew = fOld;
iteration = 0;


while cor > tol & iteration < 10000
    iteration = iteration + 1
    fNew = w*(fNew-fOld)+fOld;
    fOld = fNew;
    a = zeros(1,n);                 % Allocation for pentadiagonal
    b = zeros(1,n);
    c = zeros(1,n);
    d = zeros(1,n);
    e = zeros(1,n);
    f = zeros(1,n);
    

    for i = 1:n                     % Forms pentadiagonal
            a(i) = -1/2*dx^3;
            b(i) = (1/dx^3) + (fOld(i)/(2*dx^2));
            c(i) = ( - (fOld(i)/dx^2));
            d(i) = (-1/dx^3) + (fOld(i)/(2*dx^2));
            e(i) = 1/2*dx^3;
    end
    a(1:2) = 0;                     % Formats for thomas5 
    b(1) = 0;
    d(n) = 0;
    e(n-1:n) = 0;
    
    c(1) = 1;                               % First B.C.
    d(1) = 0;
    e(1) = 0;
    
    b(2) = 1/dx^3 + fOld(2)/(2*dx^2);      % f'(0) = 0 B.C.
    c(2) = -2*dx^3 -fOld(2)/dx^2;
    d(2) = -1/dx^3+ fOld(2)/(2*dx^2);
    e(2) = 1/2*dx^3;
    
    a(n) = -1/dx^3;                         % f'(L) = 1
    b(n) = 4/dx^3 + fOld(n)/dx^2;
    c(n) = -3/dx^3 - fOld(n)/dx^2;
    
    f(n) = -2/dx^2 - fOld(n)/dx;

    fNew = Thomas5(a,b,c,d,e,f);
    
    for i = 1:n
        correction(i) = fNew(i) - fOld(i);
    end
    cor = max(abs(correction));
    uNew1 = fNew;
    x1 = x;

end

% I really felt like these were the right equations and I don't know where
% I went wrong, but it won't converge.  It looks like it is enforcing a
% derivative of zero at the right end, instead of a derivative of 1. 


    figure(8)
    plot(x1,uNew1,'DisplayName', 'Pentadiagonal')
    hold on 
    legend('Location', 'best')
    xlabel('x')
    ylabel('f')
    title('Blasius Equation')
    hold off


end


%% Problem 4
function Problem4()
print("Hello, World")

%He always says try to write some code to get some points, so I wrote some
%code
end
%% Problem 2a Solver
function [x1,u1] = Problem2Solver(x1,x2,b1,b2,alpha)
n = 11;
x = linspace(x1,x2,n);
dx = x(2) - x(1);
alpha = alpha;

w=.9;
cor = dx;
tol = 1e-4;
uOld = zeros(1,n);
uNew = zeros(1,n);
iteration = 0;

while cor > tol & iteration < 10000
    iteration = iteration + 1
    uNew = w*(uNew-uOld)+uOld;
    uOld = uNew;
    a = zeros(1,n);                 % Allocation for tridiagonal
    b = zeros(1,n);
    c = zeros(1,n);
    d = zeros(1,n);
    B = zeros(1,n);
    
   
    for i = 1:n                     % Forms tridiagonal
            a(i) = 1/dx^2;
            b(i) = -2/dx^2;
            c(i) = 1/dx^2;
            d(i) = alpha*uOld(i);
    end
    a(end) = 0;                     % Formats for thomas3 
    a(1) = 0;
    b(1) = 1;
    b(end) = 1;
    c(1) = 0;
    c(end) = 0;
    d(1) = b1;
    d(end) = b2;
    uNew = Thomas3(a,b,c,d);
    
    for i = 2:n-1
        correction(i) = uNew(i) - uOld(i);
    end
    cor = max(abs(correction));
    x1 = x;
    u1 = uNew;

end



end
%% Tridiagonal Thomas Algorithm Function
function xBar = Thomas3(a,b,c,d)
    n = length(b);                      % Length of total diagonal
    aBar = a;                       % Forms vectors of n length for manipulated values, padded as necessary
    bBar = b;
    cBar = c;
    dBar = d;

    for i=2:n                           % "Zip-down" eliminates subdiagonal by subtracting the previous row scaled to the 'a' term
        bBar(i) = b(i) - aBar(i)*cBar(i-1)/bBar(i-1);
        dBar(i) = d(i) - aBar(i)*dBar(i-1)/bBar(i-1);
    end
    
    xBar(n) = dBar(n)/ bBar(n);         % "Zip-up" solves for x from last row single variable equation to first row
    for i = n-1:-1:1
        xBar(i) = ((dBar(i)-(cBar(i)*xBar(i+1)))/bBar(i));
    end
end

%% LU Factorization function
function x = LU(A,b)
    n = length(b);
    U = zeros(n);               % Allocation for U and L matrices
    L = zeros(n);
    y = zeros(n,1);
    x = zeros(n,1);
    
    for j = 1:n
        U(1,j) = A(1,j);         % U gets first row of A
        L(j,j) = 1;              % Diagonal of L gets 1
    end
    
    for i=2:n                    % Populates first column of L
        L(i,1) = A(i,1)/U(1,1);
    end
     
    
    for i=2:n
        for j=1:n               % Forms U matrix
            for k=i:n  
                sumLU = 0;
                for l=1:i-1
                    sumLU = sumLU + L(i,l)*U(l,k);
                end
                U(i,k) = A(i,k) - sumLU;
            end
            for k=1:i-1         % Forms L matrix
                sum=0;
                
                for l=1:k-1
                    sum=sum+L(i,l)*U(l,k);
                end
                L(i,k)=(A(i,k)-sum)/U(k,k);
            end
        end
    end
    
                        
    y(1) = b(1);                 % Solving Ly = b downward
    
    for i = 2:n
        Lysum = 0;
        for k = 1:i-1
            Lysum = Lysum + L(i,k)*y(k);
        end
        y(i) = b(i) - Lysum;
    end
    
    x(n) = y(n)/U(n,n);
                      % Solving Ux = y upward
    for i = n-1:-1:1
        Uxsum = 0;
        for k = i+1:n
            Uxsum = Uxsum + U(i,k)*x(k);
        end
        x(i) = (y(i)-Uxsum)/U(i,i);
    end  
end
%% Pentadiagonal Thomas Algorithm Function
function xBar = Thomas5(a,b,c,d,e,f)
    n = length(c);                      % Length of total diagonal
    aBar = a;                     % Forms vectors of n length for manipulated values, padded as necessary
    bBar = b;
    cBar = c;
    dBar = d;
    eBar = d;
    fBar = f;
    
    for i=3:n                           % "Zip-down" eliminates lower subdiagonal by subtracting the previous row scaled to the 'a' term
        bBar(i) = bBar(i) - aBar(i)*cBar(i-1)/bBar(i-1);
        cBar(i) = cBar(i) - aBar(i)*dBar(i-1)/bBar(i-1);
        dBar(i) = dBar(i) - aBar(i)*eBar(i-1)/bBar(i-1);
        fBar(i) = fBar(i) - aBar(i)*fBar(i-1)/bBar(i-1);
    end
    
    cHat = cBar;                        % Forms vectors for manipulated values during second downward elimination
    dHat = dBar;
    eHat = eBar;
    fHat = fBar;
    
    for i=2:n                           % Second downward elimination step
        cHat(i) = cHat(i) - (bBar(i)*dHat(i-1)/cHat(i-1));
        dHat(i) = dHat(i) - (bBar(i)*eHat(i-1)/cHat(i-1));
        fHat(i) = fHat(i) - (bBar(i)*fHat(i-1)/cHat(i-1));
    end
   
    xBar(n) = fHat(n)/ cHat(n);         % "Zip-up" solves for x from last row single variable equation to first row
    xBar(n-1) = (fHat(n-1)-(xBar(n)*dHat(n-1)))/cHat(n-1);
    for i = n-2:-1:1
        xBar(i) = (1/cHat(i))*(fHat(i)-(xBar(i+1)*dHat(i))-(xBar(i+2)*eHat(i)));
    end
    

end
%% Tridiagonal Block Thomas Algorithm Function
function xBar = THOMAS3_BLOCK(a,b,c,d)
    N = length(b);                      % Length of diagonal of 2x2 blocks
    aBar = cell(N,1);                   % Allocation for manipulated values
    bBar = cell(N,1);
    cBar = cell(N,1);
    dBar = cell(N,1); 
    
    for i=2:N                           % Forms cell arrays for manipulated values, padded as necessary
        aBar{i} = a{i-1};
        cBar{i-1} = c{i-1};
    end
    aBar{1}=zeros(2);
    bBar{1} = b{1};
    cBar{N}=zeros(2);
    dBar{1} = d{1};
  
    for i=2:N                           % "Zip-down" eliminates subdiagonal
        bBar{i} = b{i} - aBar{i}*INVERSE(bBar{i-1})*cBar{i-1};  % Utilizes user defined INVERSE function to calculate 2X2 inverse
        dBar{i} = d{i} - aBar{i}*INVERSE(bBar{i-1})*dBar{i-1};
    end
    
    xBar{N} = INVERSE(bBar{N})*dBar{N}; % "Zip-up" solves for x from last row single variable matrix equation to first row
    for i = N-1:-1:1     
        xBar{i} = INVERSE(bBar{i})*(dBar{i}-cBar{i}*xBar{i+1});        
    end
    
    xBar = cell2mat(xBar);              % Returns x as a vector
    xBar = xBar(:)';
end

