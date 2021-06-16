%Chris Hopp
%915866326
%ENG-180 Project 5: Numerical Approximations
%11/4/2020


clc
clf
clear all
Problem1a()
Problem1b()
Problem1c()
Problem2()
Problem3()
Problem4()

%% Problem 1a
function Problem1a()

% 1.a) 
f = @(x) (sin(x)).^2;
dx = [pi/10, pi/20, pi/40, pi/80];
axis = [10,20,40,80];
location=pi/3;
exact = 2*sin(location)*(cos(location));  % Exact first derivative


for i = 1:4                                        % Calculates all first derivatives
d = FirstDerivative(f,dx(i),location,'forward',1);
res_forward1(i) = abs(d - exact);
d = FirstDerivative(f,dx(i),location,'backward',1);
res_backward1(i) = abs(d-exact);
d = FirstDerivative(f,dx(i),location,'forward',2);
res_forward2(i) = abs(d - exact);
d = FirstDerivative(f,dx(i),location,'backward',2);
res_backward2(i) = abs(d-exact);
d = FirstDerivative(f,dx(i),location,'central',2);
res_central2(i) = abs(d-exact);
d = FirstDerivative(f,dx(i),location,'cubic',1);
res_cubic(i) = abs(d-exact);
end

figure(1)
loglog(axis,res_forward1,'-o','DisplayName','1st Order Forward')
hold on
loglog(axis,res_backward1,'-s','DisplayName','1st Order Backward')
loglog(axis,res_forward2,'-^','DisplayName','2nd Order Forward')
loglog(axis,res_backward2,'-*','DisplayName','2nd Order Backward')
loglog(axis,res_central2,'-x','DisplayName','2nd Order Central')
loglog(axis,res_cubic,'-d','DisplayName','Cubic')
legend()
xlabel('\pi / dx')
ylabel('Residual')
title('First Derivative')
hold off

exact = 2 * ((cos(location))^2 - (sin(location))^2);    % Exact 2nd derivative

for i = 1:4                                          % Calculates all 2nd derivatives
d = SecondDerivative(f,dx(i),location,'forward',1);
res_2forward1(i) = abs(d - exact);
d = SecondDerivative(f,dx(i),location,'backward',1);
res_2backward1(i) = abs(d-exact);
d = SecondDerivative(f,dx(i),location,'forward',2);
res_2forward2(i) = abs(d - exact);
d = SecondDerivative(f,dx(i),location,'backward',2);
res_2backward2(i) = abs(d-exact);
d = SecondDerivative(f,dx(i),location,'central',2);
res_2central2(i) = abs(d-exact);
d = SecondDerivative(f,dx(i),location,'cubic',1);
res_2cubic(i) = abs(d-exact);
end

figure(2)
loglog(axis,res_2forward1,'-o','DisplayName','1st Order Forward')
hold on
loglog(axis,res_2backward1,'-s','DisplayName','1st Order Backward')
loglog(axis,res_2forward2,'-^','DisplayName','2nd Order Forward')
loglog(axis,res_2backward2,'-*','DisplayName','2nd Order Backward')
loglog(axis,res_2central2,'-x','DisplayName','2nd Order Central')
loglog(axis,res_2cubic,'-d','DisplayName','Cubic')
legend()
xlabel('\pi / dx')
ylabel('Residual')
title('Second Derivative')
hold off

exact = -8*cos(location)*sin(location);               % Exact 3rd derivative

for i = 1:4                                          % Calculates all 3rd derivatives
d = ThirdDerivative(f,dx(i),location,'forward');
res_3forward1(i) = abs(d - exact);
d = ThirdDerivative(f,dx(i),location,'backward');
res_3backward1(i) = abs(d-exact);
d = ThirdDerivative(f,dx(i),location,'central');
res_3central(i) = abs(d-exact);
end


figure(3)
loglog(axis,res_3forward1,'-o','DisplayName','1st Order Forward')
hold on
loglog(axis,res_3backward1,'-s','DisplayName','1st Order Backward')
loglog(axis,res_3central,'-^','DisplayName','2nd Order Central')
legend()
xlabel('\pi / dx')
ylabel('Residual')
title('Third Derivative')
hold off

exact = 8*((sin(location))^2 -(cos(location))^2);    % Exact 4th derivative

for i = 1:4                                          % Calculates all 4th derivative
d = FourthDerivative(f,dx(i),location);
res_4central(i) = abs(d-exact);
end


figure(4)
loglog(axis,res_4central,'-o','DisplayName','2nd Order Central')
legend()
xlabel('\pi / dx')
ylabel('Residual')
title('Fourth Derivative')
hold off
end

%% Problem 1.b
function Problem1b()
Fxy = @(x,y) (sin(x)).^2.*(sin(y)).^2;
locationx = pi/3;
locationy = pi/3;
dx = [pi/10, pi/40];
Fxx = @(x) Fxy(x,locationy);             % Function fixed in y


for i = 1:2                             % Second partial w.r.t x
fxx(i) = SecondDerivative(Fxx,dx(i),locationx,'central',2);
fxxyy(i) = Partialxxyy(Fxy,dx(i),dx(i),locationx, locationy);
laplace(i) = Laplacian(Fxy,dx(i),dx(i),locationx,locationy);
biharmonic(i) = Biharmonic(Fxy,dx(i),dx(i),locationx,locationy);
end

solution1b = [fxx',fxxyy',laplace',biharmonic']

end

%% Problem 1.c
function Problem1c()
f = @(x) (sin(x)).^2;
order = 10;                 % Number of Taylor series terms
dx = pi/10;                 % Spacing for numerical derivatives
a = 0;                      % Centered on 0
x = linspace(-pi,pi,201);   % Input vector for TS(x)

taylor = TaylorSeries(f,a,order,x,dx);
figure(5)
plot(x,taylor,'--','DisplayName','Taylor Series')
hold on
plot(x,f(x),'DisplayName','Actual')
ylim([-10,2])
xlabel('x')
ylabel('f(x)')
title('Sin^2(x)')
legend()

end

%% Problem 2
function Problem2()
    f = @(x) (sin(x)).^2;
    dx = [pi/3, pi/6, pi/12, pi/24];
    axis = [3,6,12,24];
    a = 0;
    b = 2*pi/3;
    
    exact = (1/24)*(8*pi + 3^(3/2));
    for i = 1:4                             % Calculates integrals by each method
        I = RectangleRule(f,dx(i),a,b);
        rectangle(i) = abs(I - exact);
        I = TrapezoidRule(f,dx(i),a,b);
        trapezoid(i) = abs(I - exact);
        I = Simpson(f,dx(i),a,b);
        simpson(i) = abs(I - exact);
        [analytical, numeric] = Cubic(f,dx(i),a,b);
        cubicAnalytical(i) = abs(analytical - exact);
        cubicNumeric(i) = abs(numeric - exact);
    end
    figure(6)

    loglog(axis,rectangle,'-o','DisplayName','Rectangle Rule');
    hold on
    loglog(axis,trapezoid,'-s','DisplayName','Trapezoid Rule');
    loglog(axis,simpson,'-^','DisplayName',"Simpson's Rule");
    loglog(axis,cubicAnalytical,'-*','DisplayName','Cubic Analytical');
    loglog(axis,cubicNumeric,'-*','DisplayName','Cubic Numeric');
    xlabel('\pi/dx')
    ylabel('Residual')
    title('\int_0^{2\pi/3} sin^2(x)')
    legend()
end

%% Problem 3
function Problem3()
dx = pi/100;
a = -pi;
b = -pi/2;
c = pi/2;
d = pi;


x1 = a:dx:b; y1 = -ones(length(x1),1);   % Splits function into piecewise definition
x2 = b:dx:c; y2 = ones(length(x2),1);
x3 = c:dx:d; y3 = -ones(length(x3),1);

a0 = (1/(2*pi)) + Trap(x1,y1) + Trap(x2,y2) + Trap(x3,y3);

for j = 1:6
    n = 2^(j-1);

    for i = 1:n                         % a_n coefficients
        ya1 = @(x) -cos(i.*x);
        ya2 = @(x) cos(i.*x);
        ya3 = @(x) -cos(i.*x);
        an(i) = (1/(pi)) * (TrapezoidRule(ya1,dx,a,b) + TrapezoidRule(ya2,dx,b,c) + TrapezoidRule(ya3,dx,c,d));
    end
    for i = 1:n                         % b_n coefficients
        yb1 = @(x) -sin(i.*x);
        yb2 = @(x) sin(i.*x);
        yb3 = @(x) -sin(i.*x);
        bn(i) = (1/(pi)) * (TrapezoidRule(yb1,dx,a,b) + TrapezoidRule(yb2,dx,b,c) + TrapezoidRule(yb3,dx,c,d));
    end

    x = a:dx:d;
    y = 0;
    for i = 1:n
        y = y + an(i)*cos(i*x) + bn(i)*sin(i*x);  % Combines full series
    end
        y = a0 + y;
        figure(7)
        plot(x,y,'DisplayName', sprintf('n = %d', n))
        hold on
        xlabel('x')
        ylabel('f(x)')
        title('Fourier Series Squarewave')
        legend()
end
        hold off


dx = pi/100;
a = -pi;
b = pi;
f = @(x) x./pi;         % Sawtooth function


a0 = (1/(2*pi)) + TrapezoidRule(f,dx,a,b);

for j = 1:6
    n = 2^(j-1);

    for i = 1:n                         % a_n coefficients
        fa = @(x) (x./pi).*cos(i.*x);
        an(i) = (1/(pi)) * TrapezoidRule(fa,dx,a,b);
    end
    for i = 1:n                         % b_n coefficients
        fb = @(x) (x./pi).*sin(i.*x);
        bn(i) = (1/(pi)) * TrapezoidRule(fb,dx,a,b);
    end
    x = a:dx:b;
    y = 0;
    for i = 1:n
        y = y + an(i)*cos(i*x) + bn(i)*sin(i*x); % Combines full series
    end
        y = a0 + y;
        figure(8)
        plot(x,y,'DisplayName', sprintf('n = %d', n))
        hold on
        xlabel('x')
        ylabel('f(x)')
        title('Fourier Series Sawtooth')
        legend()
end
    hold off
        

end

%% Problem 4
function Problem4()
    f = @(x,y) ((sin(x)).^2) * ((sin(y)).^2);
    dx = [pi/3,pi/6,pi/12,pi/24];
    dy = [pi/3,pi/6,pi/12,pi/24];
    axis = [3,6,12,24];
    exact = (64*pi^2+48*sqrt(3)*pi +27) / (576);
    ax = 0;
    bx = 2*pi/3;
    ay = 0;
    by = 2*pi/3;
    
    for i = 1:4
        fx = @(x)(sin(x)).^2;
        fy = @(y)(sin(y)).^2;
        I = TrapezoidRule(fx,dx(i),ax,bx) * TrapezoidRule(fy,dx(i),ay,by);
        resTrapezoid(i) = abs(I-exact);        
    end

        
  % Linear      
        for k=1:4
            x = ax:dx(k):bx;
            nx = length(x);
            y = ay:dx(k):by;
            ny = length(y);
            sum = 0;
            for i = 1:nx-1           % Computes and sums each prism element
                for j = 1:ny-1
                    sum = sum + (1/6)*dx(k)*dy(k)*(f(x(i),y(j))+(2* f(x(i+1),y(j))) + (2* f(x(i),y(j+1))) + f(x(i+1),y(j+1)));
                end
            end
            resLinear(k) = abs(exact-sum);
        end

        
    % Bilinear
    
    for k = 1:4
    x = ax:dx(k):bx;
    nx = length(x);
    y = ay:dx(k):by;
    ny = length(y);
    integral = 0;
    
    for i = 1:nx-1
        for j = 1:ny-1       % Finds coefficients and solves 2d integral analytically  
            a1 = (f(x(i+1),y(j)) - f(x(i),y(j)))/dx(k);
            a2 = (f(x(i),y(j+1))- f(x(i),y(j))) /  dy(k);
            a3 = (1/(dx(k)*dy(k))) * ( f(x(i+1),y(j+1)) - f(x(i),y(j)) - a1*dx(k) - a2*dy(k));
            integral = integral + (dx(k)*dy(k)*f(x(i),y(j))) + (1/2)*a1*dy(k)*dx(k)^2 + (1/2)*a2*dy(k)^2*dx(k) + (1/4)*a3*dy(k)^2*dx(k)^2;
        end
    end
   
    resBilinear(k) = abs(integral- exact);
    
    end
        figure(9)
        loglog(axis,resTrapezoid,'-o','DisplayName','Trapezoid')
        hold on
        loglog(axis,resLinear,'-s','DisplayName','Linear')
        loglog(axis,resBilinear,'-^','DisplayName','Bilinear')
        xlabel('\pi / dx')
        ylabel('Residual')
        title('\int \int sin^2(x)sin^2(y)')
        legend()
    
    
end
%% First derivative function
function derivative = FirstDerivative (f,dx,value,method,order)
    n = 4;
    x = value-n*dx:dx:value+n*dx;       % Forms x-vector with enough points
    y = f(x);                           % y at each x
    n = n+1;
    
    switch method
        case 'forward'                  % Forward difference 1st and 2nd order
            switch order
                case 1
                    derivative = (y(n+1)- y(n)) / dx;
                case 2
                    derivative = ((-3/2)*y(n) + 2*y(n+1) - y(n+2)/2)/ dx;
                otherwise
                    derivative = NaN;
            end
        case 'backward'                % Backward difference 1st and 2nd order
            switch order
                case 1
                    derivative = (y(n)- y(n-1)) / dx; 
                case 2
                    derivative = ((3/2)*y(n) - 2*y(n-1) + y(n-2)/2)/ dx;
                otherwise
                    derivative = NaN;
            end
        case 'central'                          % Central difference
            derivative = (y(n+1) - y(n-1))/ (2*dx);
            
        case 'cubic'                                % Cubic spline           
            k1 = 2*(sin(value-dx))*cos(value-dx);     % Endpoint derivatives
            k3 = 2*sin(value+dx)*cos(value+dx);           
            
            A = [1,0,0;1/6, 4/6, 1/6;0,0,1];
            d = [k1; (y(n+1)- y(n-1))/(2*dx); k3];
            k = LU(A,d);
            derivative = k(2);           
        otherwise             
            derivative = NaN;
    end
                    
    

end

%% Second derivative function
function derivative = SecondDerivative (f,dx,value,method,order)
    n = 4;
    x = value-n*dx:dx:value+n*dx;       % Forms x-vector with enough points
    y = f(x);                           % y at each x
    n = n+1;
    
    switch method
        case 'forward'                  % Forward difference 1st and 2nd order
            switch order
                case 1
                    derivative = (y(n) - 2*y(n+1) + y(n+2)) / dx^2;
                case 2
                    derivative = (2*y(n) - 5*y(n+1) + 4*y(n+2) - y(n+3)) / dx^2;
                otherwise
                    derivative = NaN;
            end
        case 'backward'                % Backward difference 1st and 2nd order
            switch order
                case 1
                    derivative = (y(n) - 2*y(n-1) + y(n-2)) / dx^2;
                case 2
                    derivative = (2*y(n) - 5*y(n-1) + 4*y(n-2) - y(n-3)) / dx^2;
                otherwise
                    derivative = NaN;
            end
        case 'central'                          % Central difference
            derivative = (y(n+1) - 2*y(n)+ y(n-1))/ (dx^2);
        case 'cubic'                                    % Cubic spline
            
            k1 = 2*(sin(value-dx))*cos(value-dx);     % Endpoint derivatives
            k3 = 2*sin(value+dx)*cos(value+dx);           
            
            A = [1,0,0;1/6, 4/6, 1/6;0,0,1];
            d = [k1; (y(n+1)- y(n-1))/(2*dx); k3];
            k = LU(A,d);
                               
            derivative = (3 * ((y(n+1)-y(n))/(dx^2))) - ((k(3) + (2*k(2)))/dx);
            derivative = derivative * 2;
            
        otherwise 
            derivative = NaN;
    end
                    
    

end

%% Third derivative function
function derivative = ThirdDerivative(f,dx,value,method)
 n = 4;
    x = value-n*dx:dx:value+n*dx;       % Forms x-vector with enough points
    y = f(x);                           % y at each x
    n = n+1;
    
    switch method
        case 'forward'                  % Forward difference
             derivative = (-y(n) +3*y(n+1) - 3*y(n+2) + y(n+3)) / dx^3;

        case 'backward'                % Backward difference
             derivative =  (y(n) -3*y(n-1) + 3*y(n-2) - y(n-3)) / dx^3;

        case 'central'                          % Central difference
            derivative = ((-1/2)*y(n-2) + y(n-1)+ - y(n+1) + (1/2)*y(n+2))/ (dx^3);
        otherwise 
            derivative = NaN;
    end
                    
end

%% Fourth derivative function
function derivative = FourthDerivative(f,dx,value)
 n = 4;
    x = value-n*dx:dx:value+n*dx;       % Forms x-vector with enough points
    y = f(x);                           % y at each x
    n = n+1;
    derivative = (y(n+2) + -4*y(n+1)+6*y(n) -4*y(n-1) + y(n-2))/ (dx^4);                  
end

%% Higher order derivative function
function derivative = HigherOrder(f,dx,value,order)

% This algorithm was adapted from one found at 
% https://archive.siam.org/books/ot98/sample/OT98Chapter1.pdf
% I believe MATLAB has a similar algorithm for computing coefficients

    n = ceil(order/2);
    x = value-n*dx:dx:value+n*dx;       % Forms x-vector with enough points
    y = f(x);                           % y at each x
    

    n = length(x);                       % Solves Vandermonde
    A = ones(n,n);
    disp = (x(:)-value);    % Displacements
    for i=2:n
        A(i,:) = (disp .^ (i-1)) ./ factorial(i-1);
    end
    b = zeros(n,1);         % RHS
    b(order+1) = 1;             % Desired derivative term remains
    c = LU(A,b);            % Solve system for coefficients
    
    derivative = sum(c'.*y); % Derivative is sum of coefficient times each y

   
    
    
end

%% 2D Partial Second Derivative
function derivative = Partialxxyy(f,dx,dy,valuex,valuey)
 n = 4;
    x = valuex-n*dx:dx:valuex+n*dx;       % Forms x-vector with enough points
    y = valuey-n*dy:dy:valuey+n*dy;
 n = n+1; 
 
 dfx = zeros(9);
    for i = 1:9
        for j = 2:8
        dfx(i,j) = (f(x(j-1),y(i)) - 2*f(x(j),y(i)) + f(x(j+1),y(i)))/(dx^2);
        end
    end
    
 derivative = (dfx(n-1,n) - 2*(dfx(n,n)) + dfx(n+1,n))/dy^2;



end

%% Laplacian Operator
function laplacian = Laplacian(f,dx,dy,valuex,valuey)
fx = @(x) f(x,valuey);
fy = @(x) f(valuex,x);

partialx = SecondDerivative(fx,dx,valuex,'central',2);  % Partial 2nd wrt x
partialy = SecondDerivative(fy,dy,valuey,'central',2);  % Partial 2nd wrt y

laplacian = partialx + partialy;

end

%% Biharmonic Operator
function biharmonic = Biharmonic(f,dx,dy,valuex,valuey)
fx = @(x) f(x,valuey);
fy = @(x) f(valuex,x);

partialx = FourthDerivative(fx,dx,valuex);      % Partial 4th w.r.t. x
partialy = FourthDerivative(fy,dy,valuey);      % Partial 4th w.r.t. y
partialxy = Partialxxyy(f,dx,dy,valuex,valuey); % Partial 2nd x 2nd y cross term

biharmonic = partialx + partialy + 2 * partialxy;
end

%% Taylor Series
function series = TaylorSeries(f,a,order,x,dx)
n = length(x);
k = order;
coefficients = zeros(k,1);
coefficients(1) = f(a);
x = x - a;          % Accounts for Taylor as well as MacLaurin series
X = zeros(n,k);     % Allocation for 
for i = 2:k         % Taylor series coefficients
    coefficients(i) = (HigherOrder(f,dx,a,(i-1)))/factorial(i-1);
end

for i = 1:n
    for j = 1:k     % Forms matrix of 1, x^2, x^3 ....
        X(i,j) = x(i)^(j-1);    
    end
end

A = X * coefficients; % Matrix of coefficient times 1, x^1, x^2, ...
series = sum(A,2);    % Forms final series at each x point  
end

%% Rectangle Rule
function integral = RectangleRule(f,dx,a,b)
x = a:dx:b;         % x vector spaning integral bounds
n = length(x);
integral = 0;

for i = 1:n-1
    integral = integral + f(x(i)); % Total heights of each segment
end
integral = integral * dx;           % Area on uniform mesh
end
 
%% Trapezoid Rule
function integral = TrapezoidRule(f,dx,a,b)
x = a:dx:b;         % x vector spaning integral bounds
n = length(x);
integral = 0;

for i = 2:n-1
    integral = integral + f(x(i)); % Total heights of each segment
end
integral = (dx/2)*(f(x(1)) + f(x(n)) + 2*integral);           % Area on uniform mesh

end
%% Vector Trapezoid
function integral = Trap(x,y) 
n = length(x);
integral = 0;
dx = x(2) - x(1);
for i = 2:(n-1)
    integral = integral + y(i); % Total heights of each segment
end
    
    integral = (dx/2)*(y(1) + y(n) + 2*integral);           % Area on uniform mesh   
        
end
%% Simpson's Rule
function integral = Simpson(f,dx,a,b)
x = a:dx:b;         % x vector spaning integral bounds
n = length(x);
sum1 = 0;
sum2 = 0;

for i = 1:((n-1)/2)
    sum1 = sum1 + f(x(2*i)); 
end
for i = 1:((n-3)/2)
    sum2 = sum2 + f(x(1+2*i)); 
end
integral = (dx/3)*(f(x(1)) + f(x(n)) + 4*sum1 + 2*sum2);    % Area on uniform mesh

end

%% Cubic Integral

function [analytical, numeric] = Cubic(f,dx1,a,b)
x = a:dx1:b;
y = f(x);

n= length(y);                     % Length of input vector
k=zeros(n,1);                     % Allocation for slope vector
m=zeros(n,1);                     % Allocation for d2y/dx2
d=zeros(n,1);                     % Allocation for d3y/dx3
Sx=zeros(11,n-1);                 % Allocation for interpolated x vector
Sy= zeros(10,n-1);                % Allocation for interpolated y vector
Sy(1,1) = y(1);

k = 2.*sin(x).*cos(x);

for i = 1:n-1
    dx=(x(i+1)-x(i));               % Calculates m and d at each point
    dy=(y(i+1)-y(i));
    m(i) = ((3*dy)/(dx^2)) - ((k(i+1)+2*k(i))/dx);   
    d(i) = ((k(i+1)+k(i))/((dx)^2)) - ((2*(y(i+1)-y(i)))/((dx)^3));
    
    Sx(:,i) = linspace(x(i),x(i+1),11);     % X vector with 10 subdivisions
    
    for j=1:10                              % Forms spline points
        Sy(j,i)= y(i) + k(i)*(Sx(j,i)-x(i)) + m(i)*((Sx(j,i)-x(i))^2)+d(i)*((Sx(j,i)-x(i))^3);
    end
end

Sx(11,:) = [];                      % Eliminates duplicated points
Sx = [Sx(:);x(n)];                  % Resizes, includes (x_n,y_n) points        
Sy = [Sy(:);y(n)];

analytical = 0;                     % Analytical integral over each spline
for i = 1:n-1
    analytical = analytical + ((d(i)/4)*dx^4) + ((m(i)/3)*dx^3) + ((k(i)/2)*dx^2) + y(i)*dx; 
end

n = length(Sx);
dx = Sx(2)-Sx(1);
numeric = 0;
for i = 2:n-1                  % Numeric integral using trapezoid rule
    numeric = numeric + Sy(i); % Total heights of each segment
end
    numeric = dx/2 * (Sy(1) + Sy(n) + 2*numeric);           % Area on uniform mesh
 


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

