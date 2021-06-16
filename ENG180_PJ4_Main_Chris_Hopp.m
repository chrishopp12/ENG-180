%Chris Hopp
%915866326
%ENG-180 Project 4: Iterative Methods
%10/27/2020


clc
clf
clear all

Problem1()
Problem2()
Problem3()
Problem4a()
Problem4b()
 
%% Problem 1
function Problem1()
for j=1:2                      % Builds matrix system
n = 10^j;
a = ones(n,1); a(1) = 0;                    
b = -2*ones(n,1);
c = ones(n,1); c(n) = 0;
f = zeros(n,1);
for i=1:n
    f(i) = sin((i*pi)/(n+1));
end

                                % Solves using appropriate methods
xGuess = ones(n,1);
epsilon = 1e-10;
w = 1.5;
[iJacobi, resJacobi, corJacobi, tJacobi] = Jacobi(a,b,c,f,xGuess,epsilon);
[iGauss, resGauss, corGauss, tGauss] = Gauss(a,b,c,f,xGuess,epsilon);
[iSOR, resSOR, corSOR, tSOR] = SOR(a,b,c,f,xGuess,epsilon,w);
iJac = 1:iJacobi;
iG = 1:iGauss;
iS = 1:iSOR;

figure(j)
subplot(1,2,1)
loglog(iJac, resJacobi,'DisplayName',['Jacobi: Time is ', num2str(tJacobi),'s'])
hold on
loglog(iG, resGauss,'DisplayName',['Gauss-Seidel: Time is ', num2str(tGauss),'s'])
loglog(iS, resSOR,'DisplayName',['SOR \omega = ',num2str(w),' : Time is ', num2str(tSOR),'s'])

legend('Location','best')                              
title(['Residual n = ',num2str(n),' \epsilon = ',num2str(epsilon)])
xlabel('Number of Iterations') 
ylabel('Residual')
hold off
subplot(1,2,2)
loglog(iJac, corJacobi,'DisplayName',['Jacobi: Time is ', num2str(tJacobi),'s'])
hold on
loglog(iG, corGauss,'DisplayName',['Gauss-Seidel: Time is ', num2str(tGauss),'s'])
loglog(iS, corSOR,'DisplayName',['SOR \omega = ',num2str(w),' : Time is ', num2str(tSOR),'s'])
legend('Location','best')                              
title(['Correction n = ',num2str(n),' \epsilon = ',num2str(epsilon)])
xlabel('Number of Iterations') 
ylabel('Correction')


hold off


end

end

%% Problem 2
function Problem2()

% Solid Propellant Burning

for j = 1:2                     % Forms each case
    switch j
        case 1
            n = .5;
        otherwise
            n = 2;
    end
c1 = 1;
c2 = 1.5;
f = @(x) c1.*x.^n - c2.*x;      % Given function
f1 = @(x) c1.*x.^n;            % Function split to plot
f2 = @(x) c2.*x;
epsilon = 1e-10;
xGuess = 1.4;
xGuess2 = .2;
[xFixed,iFixed,resFixed,corFixed, tFixed] = Fixed(f,xGuess,epsilon);
[xNewton,iNewton,resNewton,corNewton, tNewton] = Newton(f,xGuess,epsilon);
[xSecant,iSecant,resSecant,corSecant, tSecant] = Secant(f,xGuess,xGuess2,epsilon);
iFixed = 1:iFixed;
iNewton = 1:iNewton;
iSecant = 1:iSecant;
xFixed
xNewton
xSecant
solution = [xFixed, xNewton, xSecant];
figure(3)
subplot(1,2,j)
fplot(f1,[-5,5])
hold on
fplot(f2,[-5,5])
xline(xGuess,'DisplayName','Initial Guess','LineStyle',':')
scatter(solution,f1(solution),'DisplayName', 'Solution')
xlabel('x')
ylabel('f(x)')
title(['Solid Propellant Burning n=',num2str(n)])
legend
hold off

figure(3+j)
subplot(1,2,1)
loglog(iFixed, resFixed,'DisplayName',['Fixed Point: Time is ', num2str(tFixed),'s'])
hold on
loglog(iNewton, resNewton,'DisplayName',['Newton Method: Time is ', num2str(tNewton),'s'])
loglog(iSecant, resSecant,'DisplayName',['Secant Method: Time is ', num2str(tSecant),'s'])

legend('Location','best')                              
title(['Residual n=',num2str(n)])
xlabel('Number of Iterations') 
ylabel('Residual')
hold off

subplot(1,2,2)
loglog(iFixed, corFixed,'DisplayName',['Fixed Point: Time is ', num2str(tFixed),'s'])
hold on
loglog(iNewton, corNewton,'DisplayName',['Newton Method: Time is ', num2str(tNewton),'s'])
loglog(iSecant, corSecant,'DisplayName',['Secant Method: Time is ', num2str(tSecant),'s'])
legend('Location','best')                              
title(['Correction n=',num2str(n)])
xlabel('Number of Iterations') 
ylabel('Correction')
hold off
end

% Semenov's Thermal Analysis

a = -0.1;
b = 0.4;
c = 2.5;    
alpha = 2;            
s = @(x) a + b.*x - c.*exp((-alpha./x));  % Given function
epsilon = 1e-10;
xGuess = .5;
xGuess2 = .2;
[xFixed,iFixed,resFixed,corFixed, tFixed] = Fixed(s,xGuess,epsilon);
[xNewton,iNewton,resNewton,corNewton, tNewton] = Newton(s,xGuess,epsilon);
[xSecant,iSecant,resSecant,corSecant, tSecant] = Secant(s,xGuess,xGuess2,epsilon);
iFixed = 1:iFixed;
iNewton = 1:iNewton;
iSecant = 1:iSecant;


xFixed
xNewton
xSecant
solution = [xFixed, xNewton, xSecant];

figure(6)
fplot(s,[0,5])
hold on
xline(xGuess,'DisplayName','Initial Guess','LineStyle',':')
scatter(solution,s(solution),'DisplayName', 'Solution')
xlabel('x')
ylabel('s(x)')
title('Semenov Thermal Analysis')
legend
hold off

figure(7)
subplot(1,2,1)
loglog(iFixed, resFixed,'DisplayName',['Fixed Point: Time is ', num2str(tFixed),'s'])
hold on
loglog(iNewton, resNewton,'DisplayName',['Newton Method: Time is ', num2str(tNewton),'s'])
loglog(iSecant, resSecant,'DisplayName',['Secant Method: Time is ', num2str(tSecant),'s'])

legend('Location','best')                              
title('Residual')
xlabel('Number of Iterations') 
ylabel('Residual')
hold off

subplot(1,2,2)
loglog(iFixed, corFixed,'DisplayName',['Fixed Point: Time is ', num2str(tFixed),'s'])
hold on
loglog(iNewton, corNewton,'DisplayName',['Newton Method: Time is ', num2str(tNewton),'s'])
loglog(iSecant, corSecant,'DisplayName',['Secant Method: Time is ', num2str(tSecant),'s'])
legend('Location','best')                              
title('Correction')
xlabel('Number of Iterations') 
ylabel('Correction')
hold off


end

%% Problem 3
function Problem3()
X = [0.1, 1.05, 1.1, 1.7];
Y = [1.3, 0, 0.97, 1.2];
Z = [0, 1.99, 1.04, 0];
Delta = [1.35, 1.76, 1.01, 1.418];

fx = cell(4,1);
fy = cell(4,1);
fz = cell(4,1);
fdelta = cell(4,1);

for i = 1:4                     % Builds system of functions
fx{i} = @(x)((x-X(i)).^2);
fy{i} = @(y)((y-Y(i)).^2);
fz{i} = @(z)((z-Z(i)).^2);
fdelta{i} = @(delta)(-((delta + Delta(i)).^2));
end
F = {fx,fy,fz,fdelta};          % Gives ability to index f_i

func = @(x,y,z,delta) (x-X).^2 + (y-Y).^2 + (z-Z).^2 - (delta + Delta).^2;
xGuess = [1,1,1,1]';
epsilon = 10e-10;

[x,iteration,residual,correction,time] = Vector(F,func,xGuess,epsilon)
iteration = 1:iteration;

figure(8)
subplot(1,2,1)
loglog(iteration, residual,'DisplayName',['Time is ', num2str(time),'s'])
hold on
legend('Location','best')                              
title('GPS Residual')
xlabel('Number of Iterations') 
ylabel('Residual')
hold off

subplot(1,2,2)
loglog(iteration, correction,'DisplayName',['Time is ', num2str(time),'s'])
hold on
legend('Location','best')                              
title('GPS Correction')
xlabel('Number of Iterations') 
ylabel('Correction')
hold off



end

 
 %% Problem 4a
function Problem4a()
tStart = tic;
f = @(z) z^4 - 1;
n = 501;
epsilon = 10e-2;
index = 0;
rx=[];ry=[];bx=[];by=[];gx=[];gy=[];yx=[];yy=[];      % Color vectors
for i= -1:(2/(n-1)): 1
    for r = 1:-(2/(n-1)):-1
        index = index + 1;
        if abs(i) ~= abs(r)
            iteration = 0;
            zNew = complex(r,i);
             zOld = zNew; 
            while ((abs(f(zOld)) > epsilon) && (iteration < 100))
                zOld = zNew;
                zNew = 3.*zOld./4 + 1./(4*zOld.^3); % old - f(old)/d(old)/dx
                iteration = iteration + 1;
                if abs(f(zNew))< epsilon                % Creates color vectors
%                     color(index) = zNew;
                    if (real(zNew) >= (1 - epsilon) & real(zNew) <= (1 + epsilon))
                        rx = [rx;r];
                        ry = [ry;i];
                    elseif (real(zNew) >= (-1 - epsilon) & real(zNew) <= (-1 + epsilon))
                        yx = [yx;r];
                        yy = [yy;i];
                    elseif (imag(zNew) >= (1 - epsilon) & imag(zNew) <= (1 + epsilon))
                        gx = [gx;r];
                        gy = [gy;i];
                    elseif (imag(zNew) >= (-1 - epsilon) & imag(zNew) <= (-1 + epsilon))
                        bx = [bx;r];
                        by = [by;i];
                         end
            
                       
                end   
            end
        end
     end
end
figure(9)
scatter(rx,ry,24,'.','r','DisplayName','z = 1')
hold on
xlabel('Re')
ylabel('Im')
title(['Basins of attraction z^4 = 1, n = ',num2str(n),' Complex Analysis'])
scatter(yx,yy,24,'.','y','DisplayName','z = -1')
scatter(gx,gy,24,'.','g','DisplayName','z = i')
scatter(bx,by,24,'.','b','DisplayName','z = -i')
legend('Location','best')
hold off

% I'm leaving all this code so you can see the nightmare that this problem
% was for me.  I originally stored each solution in a map and then used
% logical indexing to convert that map to an rgb color scheme.  Ultimately
% it let me down when the size of the plot was the size of n, not the real
% imaginary axis that I wanted.


% color
% red =  (real(color) >= (1 - epsilon) & real(color) <= (1 + epsilon))
% yellow = (real(color) >= (-1 - epsilon) & real(color) <= (-1 + epsilon))
% green = (imag(color) >= (1 - epsilon) & imag(color) <= (1 + epsilon))
% blue = (imag(color) >= (-1 - epsilon) & imag(color) <= (-1 + epsilon))
% len = length(red)
% 
% C = zeros(len,len,3);
% C(:,:,1) = red + yellow;
% C(:,:,2) = yellow + green;
% C(:,:,3) = blue;
% 
% image(C)
% figure(8)
% hold on
% map1 = [1,0,0;0,0,0];
% map2 = [1,1,0;0,0,0];
% map3 = [0,1,0;0,0,0];
% map4 = [0,0,1;0,0,0];
% 
% 
% p1 = image(red);
% p2 = image(yellow);
% p3 = image(green);
% p4 = image(blue);
% 
% 
% colormap(p1,map1)
% colormap(p2,map2)
% colormap(p3,map3)
% colormap(p4,map4)
% 
% hold off
timeComplex = toc(tStart)
 end

%% Problem 4b
function Problem4b()
tStart = tic;
J = zeros(2);
d = zeros(2,1);
residual = zeros(2,1);
xOld = zeros(2,1);
xNew = xOld;
n = 501;
epsilon = 10e-2;
index = 0;
rx=[];ry=[];bx=[];by=[];gx=[];gy=[];yx=[];yy=[];
for i= -1:(2/(n-1)): 1
    for r = 1:-(2/(n-1)):-1
        index = index + 1;
        if abs(i) ~= abs(r)
            iteration = 0;
            xNew(1) = r;
            xNew(2) = i;
            res = 2*epsilon;
            while (res > epsilon) && (iteration < 100)
                xOld = xNew;
                J(1,1) = 4*xOld(1)^3-12*xOld(1)*xOld(2)^2;      % Builds Jacobian
                J(1,2) = 4*xOld(2)^3-12*xOld(2)*xOld(1)^2;
                J(2,1) = (12*(xOld(1)^2)*(xOld(2))) - (4*xOld(2)^3);
                J(2,2) = (4*xOld(1)^3) - (12*(xOld(2)^2)*(xOld(1)));
                d(1) = xOld(1)^4 - 6*((xOld(1)^2)*(xOld(2)^2)) + xOld(2)^4 - 1;
                d(2) = 4*(xOld(1)^3)*(xOld(2)) - 4*(xOld(2)^3)*(xOld(1));
                d = -d;
                iteration = iteration + 1;
                correction=LU(J,d);                              % Solves J dx = -f
                xNew = xOld + correction;
                residual(1) = xNew(1)^4 - 6*((xNew(1)^2)*(xNew(2)^2)) + xNew(2)^4 - 1;
                residual(2) = 4*(xNew(1)^3)*(xNew(2)) - 4*(xNew(2)^3)*(xNew(1));
                res = max(abs(residual));
                if res < epsilon
                    zNew = complex(xNew(1),xNew(2));         % Creates color vectors
                    if (real(zNew) >= (1 - epsilon) & real(zNew) <= (1 + epsilon))
                        rx = [rx;r];
                        ry = [ry;i];
                    elseif (real(zNew) >= (-1 - epsilon) & real(zNew) <= (-1 + epsilon))
                        yx = [yx;r];
                        yy = [yy;i];
                    elseif (imag(zNew) >= (1 - epsilon) & imag(zNew) <= (1 + epsilon))
                        gx = [gx;r];
                        gy = [gy;i];
                    elseif (imag(zNew) >= (-1 - epsilon) & imag(zNew) <= (-1 + epsilon))
                        bx = [bx;r];
                        by = [by;i];
                    end
                end
            end
        end
    end
end
figure(10)
scatter(rx,ry,24,'.','r','DisplayName','z = 1')
hold on
xlabel('Re')
ylabel('Im')
title(['Basins of attraction z^4 = 1, n = ',num2str(n),' Real analysis'])
scatter(yx,yy,24,'.','y','DisplayName','z = -1')
scatter(gx,gy,24,'.','g','DisplayName','z = i')
scatter(bx,by,24,'.','b','DisplayName','z = -i')
legend('Location','best')
hold off
timeReal = toc(tStart)
end
%% Jacobi Method
function [iterations,residual,correction, time] = Jacobi(a,b,c,f,xGuess,epsilon)
tStart = tic;               
xNew = xGuess;
n = length(f);
res = 2*epsilon;            
iteration = 0;
while res > epsilon                         % Loop runs until convergence
    residue = zeros(n,1);
    xOld = xNew;                            % Resets xOld to next level
    
    xNew(1) = (f(1) - c(1)*xOld(2))/b(1);   % Builds xNew from xOld 
    for i=2:n-1
        xNew(i) = (f(i) - (a(i) * xOld(i-1)) - c(i)*xOld(i+1))/b(i);
    end
    xNew(n) = (f(n) - a(n)*xOld(n-1))/b(n);
    
                                           % Calculates term by term residual by Ax-f
    residue(1) = b(1)*xNew(1) + c(1)*xNew(2)-f(1);
    for i = 2:n-1
        residue(i) = abs(a(i)*xNew(i-1) + b(i)*xNew(i) + c(i)*xNew(i+1)-f(i));
    end
    residue(1) = abs(b(1)*xNew(1) + c(1)*xNew(2) - f(1));
    residue(n) = abs(a(n)*xNew(n-1) + b(n)*xNew(n) - f(n));
    res = max(residue);                     % Calculates total residual and correction 
    cor = max(abs(xNew-xOld));
    iteration = iteration + 1;              % Number of total iterations
    resPlot(iteration) = res;               % Stores residual and correction to plot
    corPlot(iteration) = cor;
end
time = toc(tStart);
iterations = iteration;
residual = resPlot;
correction = corPlot;
    
    
end

%% Gauss-Seidel Method
function [iterations,residual,correction, time] = Gauss(a,b,c,f,xGuess,epsilon)
tStart = tic;               
xNew = xGuess;
n = length(f);
res = 2*epsilon;            
iteration = 0;
while res > epsilon                         % Loop runs until convergence
    residue = zeros(n,1);
    xOld = xNew;                            % Resets xOld to next level
    
    xNew(1) = (f(1) - c(1)*xOld(2))/b(1);   % Builds xNew from xOld and xNew(i-1)
    for i=2:n-1
        xNew(i) = (f(i) - (a(i) * xNew(i-1)) - c(i)*xOld(i+1))/b(i);
    end
    xNew(n) = (f(n) - a(n)*xNew(n-1))/b(n);
    
                                           % Calculates term by term residual by Ax-f
    residue(1) = b(1)*xNew(1) + c(1)*xNew(2)-f(1);
    for i = 2:n-1
        residue(i) = abs(a(i)*xNew(i-1) + b(i)*xNew(i) + c(i)*xNew(i+1)-f(i));
    end
    residue(1) = abs(b(1)*xNew(1) + c(1)*xNew(2) - f(1));
    residue(n) = abs(a(n)*xNew(n-1) + b(n)*xNew(n) - f(n));
    res = max(residue);                   % Calculates total residual and correction 
    cor = max(abs(xNew-xOld));
    iteration = iteration + 1;              % Number of total iterations
    resPlot(iteration) = res;               % Stores residue and correction to plot
    corPlot(iteration) = cor;
end
time = toc(tStart);
iterations = iteration;
residual = resPlot;
correction = corPlot;
    
    
end

%% Successive Over-Relaxation
function [iterations,residual,correction, time] = SOR(a,b,c,f,xGuess,epsilon,w)
tStart = tic;               
xNew = xGuess;
n = length(f);
res = 2*epsilon;            
iteration = 0;
while res > epsilon                         % Loop runs until convergence
    residue = zeros(n,1);
    xOld = xNew;                            % Resets xOld to next level
    
    xNew(1) = (f(1) - c(1)*xOld(2))/b(1);   % Builds xNew from xOld and xNew(i-1)
    xNew(1) = xOld(1) + w*(xNew(1) - xOld(1));
    for i=2:n-1
        xNew(i) = (f(i) - (a(i) * xNew(i-1)) - c(i)*xOld(i+1))/b(i);
        xNew(i) = xOld(i) + w*(xNew(i) - xOld(i));
    end
    xNew(n) = (f(n) - a(n)*xNew(n-1))/b(n);
    xNew(n) = xOld(n) + w*(xNew(n) - xOld(n));
    
                                           % Calculates term by term residual by Ax-f
    residue(1) = b(1)*xNew(1) + c(1)*xNew(2)-f(1);
    for i = 2:n-1
        residue(i) = abs(a(i)*xNew(i-1) + b(i)*xNew(i) + c(i)*xNew(i+1)-f(i));
    end
    residue(1) = abs(b(1)*xNew(1) + c(1)*xNew(2) - f(1));
    residue(n) = abs(a(n)*xNew(n-1) + b(n)*xNew(n) - f(n));
    res = max(residue);                    % Calculates total residual and correction 
    cor = max(abs(xNew-xOld));
    iteration = iteration + 1;              % Number of total iterations
    resPlot(iteration) = res;               % Stores residue and correction to plot
    corPlot(iteration) = cor;
end
time = toc(tStart);
iterations = iteration;
residual = resPlot;
correction = corPlot;
    
    
end

%% Fixed Point Iteration
function [x,iteration,residual,correction,time] = Fixed(f,xGuess,epsilon)
tStart = tic;
iteration = 0;
xNew = xGuess;
while abs(f(xNew)) > epsilon                % Xnew
    xOld = xNew;
    xNew = xOld + f(xOld);
    iteration = iteration + 1;
    residual(iteration) = abs(f(xNew));
    correction(iteration) = abs(xNew-xOld);
end


x=xNew;
time = toc(tStart);

end

%% Newtons Method
function [x,iteration,residual,correction,time] = Newton(f,xGuess,epsilon)
tStart = tic;
iteration = 0;
xNew = xGuess;
while abs(f(xNew)) > epsilon
    xOld = xNew;
    dx = linspace(xOld-10,xOld+10,1000);    % This is probably too fine of resolution
    dy = f(dx);
    dydx = Derivative(dx,dy,xOld);          % Numeric derivative at xOld
    if abs(dydx) > eps
        xNew = xOld - f(xOld)/dydx;
    else
        warning('Go to hell')               % Don't divide by zero or you
        x = NaN;
        break
    end
    iteration = iteration + 1;
    residual(iteration) = abs(f(xNew));
    correction(iteration) = abs(xNew-xOld);
end
x=xNew;
time = toc(tStart);
end


%% Secant Method
function [x,iteration,residual,correction,time] = Secant(f,xGuess,xGuess2,epsilon)
tStart = tic;
iteration = 0;
xNew = xGuess;
xOld = xGuess2;
while abs(f(xNew)) > epsilon
    xOlder = xOld;
    xOld = xNew;    
    dydx = (f(xOld) - f(xOlder))/(xOld - xOlder);  % Numeric derivative at xOld/Older
    if abs(dydx) > eps
        xNew = xOld - f(xOld)/dydx;
    else
        warning('Go to hell')               % Don't divide by zero or you
        x = NaN;
        break
    end
    iteration = iteration + 1;
    residual(iteration) = abs(f(xNew));
    correction(iteration) = abs(xNew-xOld);
end
x=xNew;
time = toc(tStart);
end


%% Vector Newtons Method
function [x,iteration,residual,correction,time] = Vector(F,func,xGuess,epsilon)
tStart = tic;
iteration = 0;
xNew = xGuess;
res = 2*epsilon;
while res > epsilon
    xOld = xNew;
    J = Jacobian(F,xOld(1),xOld(2),xOld(3),xOld(4));    % Builds Jacobian(xOld)
    d = -func(xOld(1),xOld(2),xOld(3),xOld(4));         % RHS
    iteration = iteration + 1;
    correction=LU(J,d);                                 % Solves J dx = -f
    xNew = xOld + correction;

    res = max(abs(func(xNew(1),xNew(2),xNew(3),xNew(4))));
    residual(iteration) = res;
    correction(iteration) = max(abs(correction));
end
x=xNew;
time = toc(tStart);
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

%% Jacobian
function j = Jacobian(F,x,y,z,delta)
n = length(F);
J = zeros(n);
A = [x,y,z,delta];

for row = 1:n
    for column = 1:n
        x1 = linspace(A(column) -10,A(column) + 10, 100);   % Vectors for derivative
        y1 = F{column}{row}(x1);
        location = A(column);                               % Location to evaluate
        J(row, column)= Derivative(x1,y1,location);         % Forms Jacobian

    end
end
j = J;
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