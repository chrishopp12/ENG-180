%Chris Hopp
%915866326
%ENG-180 Project 3: Least Squares and Curve Fitting
%10/20/2020


clc
clf
clear all
%#ok<*MHERM>



       Problem1()                           % Calls function to complete problem 1
       Problem2()                           % Calls function to complete problem 2
       Problem3()                           % Calls function to complete problem 3

%% Function to generate f(x) and g(x) on stretch grid
function [x,f,g] = FormFunctions()
fx = @(x) exp(-x.^2);               % Given f(x)
gx = @(x) (cos(x)).^2;              % Given g(x)
kStretch = 1.1;                     % Stretch grid factor
xright = zeros(6,1);                % Allocation for stretched xvector
deltaX= zeros(6,1);                 % Allocation for dx

deltaX(1) = .2573;


xright(1) = 0;
xright(6) = pi/2;

for i = 1:4                         % Creates positive xvector
    deltaX(i) = deltaX(1)*kStretch^(i-1);
    xright(i+1) = xright(i) + deltaX(i);
end

xleft = -xright(end:-1:2);          
x = [xleft;xright];                  % Creates x vector
f = fx(x);                           % Creates f vector from f(x)
g = gx(x);                           % Creates g vector from g(x)
end
%% Problem 1
function Problem1()
fx = @(x) exp(-x.^2);               % Given f(x)
gx = @(x) (cos(x)).^2;              % Given g(x)
[x,f,g] = FormFunctions();                   % Creates functions f and g

        a = LeastSquare(x,f,1,'C');          % Least squares fits for F(x)- Cramer
        linear = @(x) a(1) + a(2).*x;

        a = LeastSquare(x,f,2,'C');
        quadratic = @(x) a(1) + a(2).*x + a(3).*(x.^2);

        a = LeastSquare(x,f,3,'C');
        cubic = @(x) a(1) + a(2).*x + a(3).*(x.^2)+a(4).*(x.^3);



        figure(1)                       % Figure output
        set(gcf, 'Position', get(0, 'Screensize'));
        subplot(2,2,1)
        hold on 
        fplot(fx,[-pi/2, pi/2],'DisplayName','Actual')
        fplot(linear,[-pi/2, pi/2],'DisplayName','Linear')
        fplot(quadratic,[-pi/2, pi/2],'m','DisplayName','Quadratic', 'LineStyle', ':', 'LineWidth', 1.5)
        fplot(cubic,[-pi/2, pi/2],'g','DisplayName','Cubic', 'LineStyle','--','LineWidth', 1.5)
        scatter(x,f,'r','Marker','^','DisplayName','Data')
        hold off
        legend                              
        title({'Least Squares Fit F(x)','Cramer Rule'})
        xlabel('x') 
        ylabel('F(x)')
    
       
        
        
        
        
        a = LeastSquare(x,g,1,'C');            % Least squares fits for G(x)- Cramer
        linear = @(x) a(1) + a(2).*x;

        a = LeastSquare(x,g,2,'C');
        quadratic = @(x) a(1) + a(2).*x + a(3).*(x.^2);

        a = LeastSquare(x,g,3,'C');
        cubic = @(x) a(1) + a(2).*x + a(3).*(x.^2)+a(4).*(x.^3);



        subplot(2,2,2)                      % Figure output
        hold on 
        fplot(gx,[-pi/2, pi/2],'DisplayName','Actual')
        fplot(linear,[-pi/2, pi/2],'DisplayName','Linear')
        fplot(quadratic,[-pi/2, pi/2],'m','DisplayName','Quadratic', 'LineStyle', ':', 'LineWidth', 1.5)
        fplot(cubic,[-pi/2, pi/2],'g','DisplayName','Cubic', 'LineStyle','--','LineWidth', 1.5)
        scatter(x,f,'r','Marker','^','DisplayName','Data')
        hold off
        legend                              
        title({'Least Squares Fit G(x)','Cramer Rule'})
        xlabel('x') 
        ylabel('G(x)')
       
       


        a = LeastSquare(x,f,1,'LU');           % Least squares fits for F(x)- LU
        linear = @(x) a(1) + a(2).*x;

        a = LeastSquare(x,f,2,'LU');
        quadratic = @(x) a(1) + a(2).*x + a(3).*(x.^2);

        a = LeastSquare(x,f,3,'LU');
        cubic = @(x) a(1) + a(2).*x + a(3).*(x.^2)+a(4).*(x.^3);



                                     % Figure output
        subplot(2,2,3)
        hold on
        fplot(fx,[-pi/2, pi/2],'DisplayName','Actual')
        fplot(linear,[-pi/2, pi/2],'DisplayName','Linear')
        fplot(quadratic,[-pi/2, pi/2],'m','DisplayName','Quadratic', 'LineStyle', ':', 'LineWidth', 1.5)
        fplot(cubic,[-pi/2, pi/2],'g','DisplayName','Cubic', 'LineStyle','--','LineWidth', 1.5)
        scatter(x,f,'r','Marker','^','DisplayName','Data')
        hold off
        legend                              
        title({'Least Squares Fit F(x)','LU Factorization'})
        xlabel('x') 
        ylabel('F(x)')

        
        
        
        
        a = LeastSquare(x,g,1,'LU');          % Least squares fits for G(x)- LU
        linear = @(x) a(1) + a(2).*x;

        a = LeastSquare(x,g,2,'LU');
        quadratic = @(x) a(1) + a(2).*x + a(3).*(x.^2);

        a = LeastSquare(x,g,3,'LU');
        cubic = @(x) a(1) + a(2).*x + a(3).*(x.^2)+a(4).*(x.^3);



        subplot(2,2,4)                      % Figure output
        hold on
        fplot(gx,[-pi/2, pi/2],'DisplayName','Actual')
        fplot(linear,[-pi/2, pi/2],'DisplayName','Linear')
        fplot(quadratic,[-pi/2, pi/2],'m','DisplayName','Quadratic', 'LineStyle', ':', 'LineWidth', 1.5)
        fplot(cubic,[-pi/2, pi/2],'g','DisplayName','Cubic', 'LineStyle','--','LineWidth', 1.5)
        scatter(x,f,'r','Marker','^','DisplayName','Data')
        hold off
        legend                              
        title({'Least Squares Fit G(x)','LU Factorization'})
        xlabel('x') 
        ylabel('G(x)')
        print('-depsc', 'LeastSquares')
      

end
%% Problem 2
function Problem2()
for k=1:3                       % Loop for grid spacing
    switch k
        case 1
            n=13;
        case 2
            n=17;
        otherwise
            n=21;
    end


x = linspace(0,2*pi,n);                     % Function to generate data
f = sin(x) + .25*sin(8.5*x);
fx = @(x)sin(x);
figure(2+k)                                 % Plotting
hold on
scatter(x,f,'DisplayName','Data')
fplot(fx,[0,2*pi],'DisplayName','Actual')

for i = 1:3                                 % Loop for fitting smoothing parameter
alpha = i*.333333333333;

a = (alpha-1)*ones(n-1,1);                  % Tridiagonal system to solve smooth fhat
b = (2-alpha)*ones(n,1);
c = a;
d = alpha*f;

fhat = THOMAS3(a,b,c,d);
fhat(1) = f(1);
fhat(n) = f(n);

m = zeros(n,1);                             % m(1) m(n) from fhat 
deltaX = x(2)-x(1);
m(1)=(1/deltaX^2)*(fhat(1)-2*fhat(2)+fhat(3));
m(n)=(1/deltaX^2)*(fhat(n)-2*fhat(n-1)+fhat(n-2));

a = ((1-alpha)/deltaX^3)* ones(n-2,1);      % Pentadiagonal system for S points
b = ((4*(alpha - 1))/deltaX^3) * ones(n-1,1);
c = (((6*(1 - alpha))/deltaX^3) + alpha) * ones(n,1);
d = b;
e = a;
h = alpha*f;

S = THOMAS5(a,b,c,d,e,h);
S(1) = f(1);
S(n) = f(n);

[Sx,Sy] = NaturalCubicSpline(x,S,m(1),m(n));  % Builds spline from S,m(1), m(n)



plot(Sx,Sy,'DisplayName',['\alpha = ',num2str(i),'/3'])
xlabel('x') 
ylabel('f(x)')
legend                                  %If line for spline data does not show properly in legend, type opengl software in command window
title(['Cubic Spline w/ Smoothing N = ',num2str(n)])
print('-depsc', ['Smoothing',num2str(k)])
end
snapnow()
hold off






end

end

%% Problem 3
function Problem3()
A = [11,0,-10,4;0,7,-29,1;-10,-29,4,0;4,1,0,1];
b = [1,0,2,0]';
C = [1,0,1,0;0,1,0,1];
d = [0,0]';
x = zeros(100,4);



epsilon = linspace(0,10000);                    % Vector of epsilon values
    for j = 1:100
        Abar = A'*A + epsilon(j)*(C'*C);        % Creates (A'A+eC'C) matrix
        dbar = A'*b + epsilon(j)*(C'*d);        % Creates (A'b+eC'd) matrix
        x(j,:) = Cramer(Abar, dbar);            % solves for x
    end
        figure(6)
        set(gcf, 'Position', get(0, 'Screensize'));
    for i = 1:4                                 % Plots x1,2,3,4 vs epsilon, semilog

        subplot(1,4,i)
        semilogx(epsilon,x(:,i))
        xlabel('\epsilon') 
        ylabel(['x',num2str(i)])
        title(['Weighted Least Squares x',num2str(i)])
        print('-depsc', ['Constraint',num2str(i)])

    end
    snapnow()
    
    Abar = A'*A + (C'*C);        % Creates (A'A+eC'C) matrix with e = 1
    dbar = A'*b + (C'*d);        % Creates (A'b+eC'd) matrix with e = 1
    xWLS = [LU(Abar, dbar)]';   % solves for x
    

    lagrangeLHS = [[2*A'*A], [C'];[C],[0,0;0,0]];        % Lagrange multipliers method
    lagrangeRHS = [[2*A'*b];[d]];
    lagrangeSol = LU(lagrangeLHS, lagrangeRHS);
    xlagrange = [lagrangeSol(1:4)]';
    
    
        
    N = null(C);                                       % Null space method 
    x0 = N(:,1);
    Ahat = A*N;
    bhat = b - A*x0;
    
    A1 = Ahat'*Ahat;
    A2 = Ahat'*bhat;
    
    z = LU(A1,A2);
    xnull = x0 +N*z;
    xnull = xnull';
    
    
    classA = [21, -4;29,6;-14,-29;4,0];     % From solving constraint x1=-x3, x2=-x4
    AL = classA'*classA;
    AR = classA'*b;
    
    xClass = LU(AL,AR);
    xClass = [xClass',-xClass'];


    tableMat = [xWLS; xlagrange; xnull; xClass];           %Output data table
    x1 = tableMat(:,1);
    x2 = tableMat(:,2);
    x3 = tableMat(:,3);
    x4 = tableMat(:,4);
    T = table(x1,x2,x3,x4, 'RowNames', {'Weighted Least Squares', 'Lagrange Multipliers', 'Nullspace', 'In-class Method'})
 

    
end
%% Function to perform least squares fit
function a = LeastSquare(x,y,n,method)
n = n+1;                    % Takes order of system
A = zeros(n);               % Allocation for A, b matrices
b = zeros(n,1);

for i=1:n                   % Forms minimization matrices based on order 
    for j=1:n
        A(i,j)=sum(x.^((j-1)+(i-1)));
    end
    b(i) = sum(y.*x.^(i-1));
end

switch method
    case 'C'
        a = Cramer(A,b);                % Solves via Cramer rule
    case 'LU'
        a = LU(A,b);                    % Solves via LU factorization              
    otherwise
        a = A\b;                        % Solves via matlab function
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
%% Cramer rule function
function x = Cramer(A,d)
n = length(A);
B = zeros(n,n,n);                    % Stores modified matrices
x = zeros(n,1);
detA = Determinant(A);
for i=1:n                   
   B(:,:,i)=A;              
   B(:,i,i)=d;                      % Replaces column i with d
   x(i) = Determinant(B(:,:,i));    % Calculates determinant of ith element
end

x = x/detA;                         % Returns solution vector


end
%% Determinant function
function detA =  Determinant(A)

if length(A) == 2
detA = A(1,1)*A(2,2) - A(1,2)*A(2,1); % Calculates 2x2 determinant
else                                  % Recursion call for size > 2x2     
    for i=1:length(A)               % Recursion statement
        new = A;
        new(1,:) = [];
        new(:,i) = [];
        z(i) = A(1,i)*Determinant(new);        
    end

    for i = 2:2:length(z)
        z(i) = -1*z(i);
    end
    detA = sum(z);       % Returns determinant for size>2
end
end

%% Tridiagonal Thomas Algorithm Function
function xBar = THOMAS3(a,b,c,d)
    n = length(b);                      % Length of total diagonal
    aBar = [0;a];                       % Forms vectors of n length for manipulated values, padded as necessary
    bBar = b;
    cBar = [c;0];
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

%% Pentadiagonal Thomas Algorithm Function
function xBar = THOMAS5(a,b,c,d,e,f)
    n = length(c);                      % Length of total diagonal
    aBar = [0;0;a];                     % Forms vectors of n length for manipulated values, padded as necessary
    bBar = [0;b];
    cBar = c;
    dBar = [d;0];
    eBar = [e;0;0];
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

%% Natural cubic spline function
function [Sx, Sy] = NaturalCubicSpline(x,y,m1,mn)
n= length(y);                     % Length of input vector
m=zeros(n,1);                     % Allocation for d2y/dx2, m(1)=m(n)=0
Sx=zeros(11,n-1);                 % Allocation for interpolated x vector
Sy= zeros(10,n-1);                % Allocation for interpolated y vector
Sy(1,1) = y(1);
m(1) = m1;
m(n) = mn;

a=(1/6)*ones(n-3,1);              % Forms tridiagonal system to solve for m vector
b=(2/3)*ones(n-2,1);
c=a;
f=zeros(n-2,1);
for i=2:n-1
    f(i-1)=(y(i+1)-(2*(y(i)))+(y(i-1)))/((x(i+1)-x(i))^2);
end
mVec = THOMAS3(a,b,c,f);         % Forms m vector
m(2:n-1)=mVec;

for i = 1:n-1
    Sx(:,i) = linspace(x(i),x(i+1),11);     % X vector with 10 subdivisions
    for j=1:10                              % Forms spline points
        dx = x(i+1) - x(i);
        a0 = m(i)/(6*dx);
        a1 = m(i+1)/(6*dx);
        a2 = (y(i)/dx)-((m(i)*dx)/6);
        a3 = (y(i+1)/dx)-((m(i+1)*dx)/6);
        Sy(j,i)= (a0*((x(i+1)-Sx(j,i))^3)) + (a1*((Sx(j,i))-x(i))^3)...
            + (a2*(x(i+1)-Sx(j,i))) + (a3*((Sx(j,i))-x(i)));
    end
end

Sx(11,:) = [];                      % Eliminates duplicated points
Sx = [Sx(:);x(n)];                  % Resizes, includes (x_n,y_n) points        
Sy = [Sy(:);y(n)];
end