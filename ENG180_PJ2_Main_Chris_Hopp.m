%Chris Hopp
%915866326
%ENG-180 Project 2: Splines
%10/13/2020


clc
clf
clear all


   Problem1()                           % Calls function to complete problem 1
   Problem2()                           % Calls function to complete problem 2
   Problem3()                           % Calls function to complete problem 3
   Problem4()                           % Calls function to complete problem 4
   Problem5()                           % Calls function to complete problem 5

%% Function to generate f(x) on stretch grid
function [x,y] = FormFunctionStretch(caseNumber)
f = @(x) cos(x).*sin(x);            % Given f(x)
kStretch = 1.1;                     % Stretch grid factor

switch caseNumber                   % Creates each test case
    case 1
        n = 10;
        deltaX = zeros(n,1);
        deltaX(1) = 0.4627;
    case 2
        n = 20;
        deltaX = zeros(n,1);
        deltaX(1) = 0.1228;
    case 3
        n = 40;
        deltaX = zeros(n,1);
        deltaX(1) = 0.0157;
end

x = zeros(n,1);                     % Allocation for stretched xvector
x(1) = 0;
x(n) = 2*pi;

for i = 1:n-2                       % Creates xvector
    deltaX(i) = deltaX(1)*kStretch^(i-1);
    x(i+1) = x(i) + deltaX(i);
end
y = f(x);                           % Creates yvector from f(x)


end

%% Function to generate f(x) on uniform grid
function [x,y] = FormFunctionUniform(caseNumber)
f = @(x) cos(x).*sin(x);            % Given f(x)

switch caseNumber                   % Creates each test case
    case 1
        n = 10;
    case 2
        n = 20;
    case 3
        n = 40;
end

x = linspace(0,2*pi,n);              % Creates xvector
y = f(x);                           % Creates yvector from f(x)
end


%% Function to produce Problem 1 solution
function Problem1()
derivative=zeros(3);                % Allocation for solutions
exactDerivative = zeros(3,1);
integral=zeros(3,1);

    for i=1:3                           % Loop to go through each test case

        [x,y] = FormFunctionStretch(i);        % Generates f(x) to be interpolated
        [Sx,Sy] = LinearSpline(x,y);    % Performs linear interpolation

        integral(i) = Integrate(Sx,Sy); % Numerical integral of splines

        for j=1:3                       % Derivative at Theta1,2,3
            derivative(i,j) = Derivative(Sx,Sy,j);
            exactDerivative(j) = ((cos(j))^2)-((sin(j))^2);
        end

        figure(1)                       % Figure output
        subplot(2,3,i)
  
        plot(Sx,Sy,'DisplayName','Spline')
        hold on
        scatter(x,y,'r','DisplayName','Data')
        hold off
        legend                                    %If line for spline data does not show properly in legend, type opengl software in command window
        title(['Linear Spline Case ' num2str(i)])
        xlabel('theta') 
        ylabel('f(theta)')
        print('-depsc', 'LinearSplineFig')
    end
nMat = [10;20;40]; results = [nMat,integral,derivative];  exact = [inf,0,exactDerivative']; % Solutions
tab = [results; exact];
array2table(tab, 'VariableNames', {'n','Integral','Derivative Theta = 1','Derivative Theta = 2','Derivative Theta = 3'},'RowNames', {'Case 1','Case 2','Case 3','Actual'})
end


%% Function to produce Problem 2 solution
function Problem2()
derivative=zeros(3);                % Allocation for solutions
exactDerivative = zeros(3,1);
integral=zeros(3,1);

    for i=1:3                           % Loop to go through each test case

        [x,y] = FormFunctionStretch(i);        % Generates f(x) to be interpolated
        [Sx,Sy] = QuadraticSpline(x,y);    % Performs linear interpolation

        integral(i) = Integrate(Sx,Sy); % Numerical integral of splines

        for j=1:3                       % Derivative at Theta1,2,3
            derivative(i,j) = Derivative(Sx,Sy,j);
            exactDerivative(j) = ((cos(j))^2)-((sin(j))^2);
        end

        figure(2)                       % Figure output
        subplot(2,3,i)
        plot(Sx,Sy,'DisplayName','Spline')
        hold on
        scatter(x,y,'r','DisplayName','Data')
        hold off
        legend                                    %If line for spline data does not show properly in legend, type opengl software in command window
        title(['Quadratic Spline Case ' num2str(i)])
        xlabel('theta') 
        ylabel('f(theta)')
        print('-depsc', 'QuadraticSplineFig')

    end
nMat = [10;20;40]; results = [nMat,integral,derivative];  exact = [inf,0,exactDerivative']; % Solutions
tab = [results; exact];
array2table(tab, 'VariableNames', {'n','Integral','Derivative Theta = 1','Derivative Theta = 2','Derivative Theta = 3'},'RowNames', {'Case 1','Case 2','Case 3','Actual'})
end


%% Function to produce Problem 3 solution
function Problem3()
derivative=zeros(3);                % Allocation for solutions
exactDerivative = zeros(3,1);
integral=zeros(3,1);

    for i=1:3                           % Loop to go through each test case

        [x,y] = FormFunctionUniform(i);        % Generates f(x) to be interpolated
        [Sx,Sy] = ClampedCubicSpline(x,y);    % Performs linear interpolation

        integral(i) = Integrate(Sx,Sy); % Numerical integral of splines

        for j=1:3                       % Derivative at Theta1,2,3
            derivative(i,j) = Derivative(Sx,Sy,j);
            exactDerivative(j) = ((cos(j))^2)-((sin(j))^2);
        end

        figure(3)                       % Figure output
        subplot(2,3,i)
        plot(Sx,Sy,'DisplayName','Spline')
        hold on
        scatter(x,y,'r','DisplayName','Data')
        hold off
        legend                                    %If line for spline data does not show properly in legend, type opengl software in command window
        title(['Clamped Cubic Case ' num2str(i)])
        xlabel('theta') 
        ylabel('f(theta)')
        print('-depsc', 'ClampedCubicFig')

    end
nMat = [10;20;40]; results = [nMat,integral,derivative];  exact = [inf,0,exactDerivative']; % Solutions
tab = [results; exact];
array2table(tab, 'VariableNames', {'n','Integral','Derivative Theta = 1','Derivative Theta = 2','Derivative Theta = 3'},'RowNames', {'Case 1','Case 2','Case 3','Actual'})
end


%% Function to produce Problem 4 solution
function Problem4()
derivative=zeros(3);                % Allocation for solutions
exactDerivative = zeros(3,1);
integral=zeros(3,1);

    for i=1:3                           % Loop to go through each test case

        [x,y] = FormFunctionUniform(i);        % Generates f(x) to be interpolated
        [Sx,Sy] = NaturalCubicSpline(x,y);    % Performs linear interpolation

        integral(i) = Integrate(Sx,Sy); % Numerical integral of splines

        for j=1:3                       % Derivative at Theta1,2,3
            derivative(i,j) = Derivative(Sx,Sy,j);
            exactDerivative(j) = ((cos(j))^2)-((sin(j))^2);
        end

        figure(4)                       % Figure output
        subplot(2,3,i)
        plot(Sx,Sy,'DisplayName','Spline')
        hold on
        scatter(x,y,'r','DisplayName','Data')
        hold off
        legend                                    %If line for spline data does not show properly in legend, type opengl software in command window
        title(['Natural Cubic Case ' num2str(i)])
        xlabel('theta') 
        ylabel('f(theta)')
        print('-depsc', 'NaturalCubicFig')

    end
nMat = [10;20;40]; results = [nMat,integral,derivative];  exact = [inf,0,exactDerivative']; % Solutions
tab = [results; exact];
array2table(tab, 'VariableNames', {'n','Integral','Derivative Theta = 1','Derivative Theta = 2','Derivative Theta = 3'},'RowNames', {'Case 1','Case 2','Case 3','Actual'})
end

%% Function to produce Problem 5 solution
function Problem5()
f = @(x,y) cos(pi*x).*cos(pi*y);    % Given f(x,y)
v=0;                                % Used to index output figures
for k=1:3
    switch k
        case 1
            n=5;
        case 2
            n=9;
        otherwise
            n=17;
    end

x = linspace(0,2,n);                % Creates xvector
y = linspace(0,2,n);                % Creates yvector
z = zeros(n);
xMat = zeros(n);
yMat = zeros(n);
for i=1:n                           % Creates grid.... there has to be an easier way
    for j=1:n
        z(i,j) = f(x(i),y(j));
    end
    xMat(:,i) = x;
    yMat(i,:) = y;
end
z1=z(:);                            % Scatter points
x1 = xMat(:);
y1 = yMat(:);
SyMat = zeros(n,(10*(n-1)+1));
for i=1:n
    [Sx,Sy] = ClampedCubicSpline(x,z(:,i));
    SyMat(i,:) = Sy;
end
SzMat = zeros(length(Sx));
for i = 1:length(Sx)
    [Sy,Sz] = ClampedCubicSpline(y,SyMat(:,i));
    SzMat(:,i) = Sz;
end

        v = v+1;
        figure(k+3+v) 
        surf(Sx,Sy,SzMat,'DisplayName','Spline')
        hold on
        scatter3(x1,y1,z1,'r','filled','DisplayName','Data')
        hold off
        title(['2-D Cubic Spline Case ' num2str(k)])
        xlabel('x') 
        ylabel('y')
        zlabel('f(x,y)')
        legend                                    %If line for spline data does not show properly in legend, type opengl software in command window
        colorbar
        zlim([-1 1])
        print('-depsc', ['2DCubic',num2str(k)])
        snapnow
        
        figure(k+4+v)
        surf(Sx,Sy,SzMat,'DisplayName','Spline')
        hold on
        scatter3(x1,y1,z1,'r','filled','DisplayName','Data')
        hold off
        title(['2-D Cubic Spline Case ' num2str(k)])       
        xlabel('x') 
        ylabel('y')
        zlabel('f(x,y)')
        
        legend                                    %If line for spline data does not show properly in legend, type opengl software in command window
        colorbar
        zlim([-1 1])
        view(2)
        print('-depsc', ['2DCubicTop',num2str(k)])
        snapnow


        
end
        figure(11)              % Displays figures for actual function
        fsurf(f,[0 2])
        title('2-D Cubic Actual')
        xlabel('x') 
        ylabel('y')
        zlabel('f(x,y)')
        legend                                    %If line for spline data does not show properly in legend, type opengl software in command window
        colorbar
        zlim([-1 1])
        print('-depsc', '2DCubicActual')
        snapnow
        
        figure(12)
        fsurf(f,[0 2])
        title('2-D Cubic Actual')
        xlabel('x') 
        ylabel('y')
        zlabel('f(x,y)')
        legend                                    %If line for spline data does not show properly in legend, type opengl software in command window
        colorbar
        zlim([-1 1])
        view(2)
        print('-depsc', '2DCubicActualTop')
        snapnow
end

%% Linear spline function
function [Sx, Sy] = LinearSpline(x,y)
n= length(y);                       % Length of input vector
k=zeros(n-1,1);                     % Allocation for slope vector
Sx=zeros(11,n-1);                   % Allocation for interpolated x vector
Sy= zeros(10,n-1);                  % Allocation for interpolated y vector
Sy(1,1) = y(1);

for i = 1:n-1
    k(i)=(y(i+1) - y(i))/(x(i+1) - x(i));   % Calculates slope
    Sx(:,i) = linspace(x(i),x(i+1),11);     % X vector with 10 subdivisions
    for j=1:10                              % Forms spline points
        Sy(j,i)= y(i)+k(i)*(Sx(j,i)-x(i));
    end
end

Sx(11,:) = [];                      % Eliminates duplicated points
Sx = [Sx(:);x(n)];                  % Resizes, includes (x_n,y_n) points        
Sy = [Sy(:);y(n)];


end




%% Quadratic spline function
function [Sx, Sy] = QuadraticSplineLeft(x,y)

n= length(y);                     % Length of input vector
k=zeros(n-1,1);                   % Allocation for slope vector
m=zeros(n-1,1);                   % Allocation for d2y/dx2
Sx=zeros(11,n-1);                % Allocation for interpolated x vector
Sy= zeros(10,n-1);               % Allocation for interpolated y vector
Sy(1,1) = y(1);

deltaX2 = x(2) - x(1);            % Calculates k1 and m1 from boundary conditions
deltaX3 = x(3) - x(1);
km = Cramer([deltaX2, deltaX2^2; deltaX3, deltaX3^2],[y(2)-y(1);y(3)-y(1)]);
k(1)=km(1);


for i=1:n-1                       % Calculates m and k values  
    m(i) = (y(i+1)-y(i)-k(i)*(x(i+1)-x(i)))/((x(i+1)-x(i))^2);
    k(i+1) = k(i) + 2*m(i)*(x(i+1)-x(i));
end 

for i = 1:n-1
    Sx(:,i) = linspace(x(i),x(i+1),11);    % X vector with 10 subdivisions
    for j=1:10                              % Forms spline points - left bias
        Sy(j,i)= y(i) + k(i)*(Sx(j,i)-x(i)) + m(i)*((Sx(j,i)-x(i))^2);
    end
end

Sx(11,:) = [];                       % Eliminates duplicated points
Sx = [Sx(:);x(n)];                  % Resizes, includes (x_n,y_n) points        
Sy = [Sy(:);y(n)];
end


function [Sx,Sy]= QuadraticSpline(x,y)

[SxLeft, SyLeft] = QuadraticSplineLeft(x,y);    % Calculates spline with left bias

x=flipud(x);                          % Flips x,y to utilize same code for right bias
y=flipud(y);

[~, SyRight] = QuadraticSplineLeft(x,y);    % Solves right bias
SyRight = flipud(SyRight);                  % Restores order of axis

Sx = SxLeft;                                % Outputs averages of R/L
Sy = (SyRight + SyLeft)/2;

end


%% Clamped cubic spline function
function [Sx, Sy] = ClampedCubicSpline(x,y)
n= length(y);                     % Length of input vector
k=zeros(n,1);                     % Allocation for slope vector
m=zeros(n,1);                     % Allocation for d2y/dx2
d=zeros(n,1);                     % Allocation for d3y/dx3
Sx=zeros(11,n-1);                 % Allocation for interpolated x vector
Sy= zeros(10,n-1);                % Allocation for interpolated y vector
Sy(1,1) = y(1);

deltaX2 = x(2) - x(1);            % Calculates k1 from boundary conditions
deltaX3 = x(3) - x(1);
deltaX4 = x(4) - x(1);
kmd = Cramer([deltaX2, deltaX2^2, deltaX2^3; deltaX3, deltaX3^2, deltaX3^3;...
    deltaX4, deltaX4^2, deltaX4^3],[y(2)-y(1);y(3)-y(1);y(4)-y(1)]);
k(1)=kmd(1);

deltaN1 = x(n-1) - x(n);           % Calculates kn from boundary conditions
deltaN2 = x(n-2) - x(n);
deltaN3 = x(n-3) - x(n);
kn = Cramer([deltaN1, deltaN1^2, deltaN1^3; deltaN2, deltaN2^2, deltaN2^3;...
    deltaN3, deltaN3^2, deltaN3^3],[y(n-1)-y(n);y(n-2)-y(n);y(n-3)-y(n)]);
k(n)=kn(1);

a=(1/6)*ones(n-3,1);                % Forms tridiagonal system to solve for k vector
b=(2/3)*ones(n-2,1);
c=a;
f=zeros(n-2,1);
for i=2:n-1
    f(i-1)=(y(i+1)-(y(i-1)))/(2*(x(i+1)-x(i)));
end
kVec = THOMAS3(a,b,c,f);            % Forms k vector
k(2:n-1)=kVec;

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
end


%% Natural cubic spline function
function [Sx, Sy] = NaturalCubicSpline(x,y)
n= length(y);                     % Length of input vector
m=zeros(n,1);                     % Allocation for d2y/dx2, m(1)=m(n)=0
Sx=zeros(11,n-1);                 % Allocation for interpolated x vector
Sy= zeros(10,n-1);                % Allocation for interpolated y vector
Sy(1,1) = y(1);

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

%% Numeric derivative function
function m = Derivative(x,y,value)
n = length(y);
location = find(x>=value,1,'first');

switch location
    case 1                              % Forward difference
        m = (y(2)-y(1))/(x(2)-x(1));
    case n                              % Backward difference
        m = (y(n)-y(n-1))/(x(n)-x(n-1));
    otherwise                           % Central difference
        m = (y(location+1)- y(location-1))/(2*(x(location)- x(location-1)));
end
end

%% Numeric integral function
function area = Integrate(x,y)
n = length(y);
area = 0;
for i=1:n-1                             % Trapezoid rule for non-uniform x
    area = area + (x(i+1)-x(i))*(y(i)+y(i+1))/2;
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
    positive = 0;                     % Initialization of positive elements                   
    negative = 0;                     % Initialization of negative elements
    for i=1:2:length(A)               % Recursion statement for odd
        new = A;
        new(1,:) = [];
        new(:,i) = [];
        positive = positive + A(1,i)*Determinant(new);        
    end
    
    for i=2:2:length(A)              % Recursion statement for even
        new = A;
        new(1,:) = [];
        new(:,i) = [];
        negative = negative - A(1,i)*Determinant(new);        
    end
    detA = positive+negative;       % Returns determinant for size>2
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

