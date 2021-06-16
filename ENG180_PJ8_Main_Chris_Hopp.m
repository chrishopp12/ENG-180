%Chris Hopp
%915866326
%ENG-180 Project 8: Partial Differential Equations
%11/29/2020


clc
clf
clear all




Problem1a()
Problem1b()
Problem2a()
Problem2b()



%% Problem 1a
function Problem1a()
x = linspace(0,pi,101);
dx = x(2) - x(1);
dt = .99*((dx^2)/2);                % Stability restriction
t = [0:dt:pi];
k = 1;
nx = length(x);
nt = length(t);


% Forward Euler
u = zeros(nx,nt);   % u(x,t)
u(1,:) = 0;
u(nx,:) = 0;
u(1:nx,1) = sin(x(1:nx));

for n = 1:nt-1                          % March in t
    for i = 2:nx-1                      % March in x 
        u(i,n+1) = ((k*dt/(dx^2))* (u(i-1,n) - 2*u(i,n) + u(i+1,n))) + u(i,n);
    end
end

figure(1)
subplot(2,2,1)
imagesc(t,x,u)
set(gca,'YDir','normal')
colorbar
xlabel('t')
ylabel('x')
title('Heat Equation: Forward Euler')


% Backward Euler
u = zeros(nx,nt);   % u(x,t)
u(1,:) = 0;
u(nx,:) = 0;
u(1:nx,1) = sin(x(1:nx));

for n = 1:nt-1                          % March in t
a = ones(nx-1,1)*(-k*dt/dx^2);
b = (ones(nx,1)*(2*k*dt/dx^2))+1;
c = a;
d = u(:,n);

sol = Thomas3(a,b,c,d);
u(1:nx,n+1) = sol(1:nx);

end
subplot(2,2,2)
imagesc(t,x,u)
set(gca,'YDir','normal')
colorbar
xlabel('t')
ylabel('x')
title('Heat Equation: Backward Euler')



% Trapezoid Rule
u = zeros(nx,nt);   % u(x,t)
u(1,:) = 0;
u(nx,:) = 0;
u(1:nx,1) = sin(x(1:nx));

for n = 1:nt-1                          % March in t
a = ones(nx-3,1)*(-k*dt/(2*dx^2));
b = (ones(nx-2,1)*(k*dt/(dx^2)))+1;
c = ones(nx-3,1)*(-k*dt/(2*dx^2));
    for i = 2:nx-1                      % March in x
        d(i-1) = u(i,n) + ((k*dt/(2*dx^2)) *(u(i-1,n)-2*(u(i,n))+u(i+1,n)));
    end

sol = Thomas3(a,b,c,d);
u(2:nx-1,n+1) = sol(1:end);

end
subplot(2,2,3)
imagesc(t,x,u)
set(gca,'YDir','normal')
colorbar
xlabel('t')
ylabel('x')
title('Heat Equation: Trapezoid Rule')



% Runge-Kutta 4
u = zeros(nx,nt);   % u(x,t)
u(1,:) = 0;
u(nx,:) = 0;
u(1:nx,1) = sin(x(1:nx));
                                    %Vector function of all u(j) derivatives
f = @(u1)[0;(u1(1:end-2)-2.*u1(2:end-1)+u1(3:end))/dx^2;0];

for j = 1:nt-1
    u(:,j+1) = RK4(f,u(:,j),dt);
end

subplot(2,2,4)
imagesc(t,x,u)
set(gca,'YDir','normal')
colorbar
xlabel('t')
ylabel('x')
title('Heat Equation: Runge-Kutta 4')


end

%% Problem 1b
function Problem1b()      % It really didn't want to cooperate to make plots
for i = 1:3
    Re = [100,1000,10000];
Problem1bCalcs(Re(i),i)
end
end
function Problem1bCalcs(Re,q)

dx = 0.01;
dy = 0.00001;
x = [-2*dx:dx:4];
y = [0:dy:10/sqrt(Re)];
nx = length(x);
ny = length(y);


u = zeros(ny,nx);
v = zeros(ny,nx);

u(:,1) = 1;
u(ny,:) = 1;


tol = 1e-2;



for j = 2:nx
    cor = 1;
    uNew = u(:,j);
    vNew = v(2:ny,j);
    uOld = ones(size(uNew));
    vOld = zeros(size(vNew));
    iteration = 0;
while cor>tol
    iteration = iteration+1;
    w=1.1;                                  % SOR
    uOld = w*uNew + (1-w)*uOld;
    vOld = w*vNew + (1-w)*vOld;
    
    if x(j) <= 0                            % Locations before plate

        a = ones(ny-1,1)*(1/(Re*dy^2));     % Length formatted for Thomas3      
        a(end) = 0;                         % Last row [0,... 0,1][u_ny] = 1
        b = (-2/(Re*dy^2)) - uOld/dx;
        b(end) = 1;
        c = ones(ny-1,1)*(1/(Re*dy^2));
        c(1) = 2*c(1);                      % u1 using ficticous point
        for k = 1:ny
            d(k) = (-uOld(k)*(u(k,j-1)))/dx;
        end
        d(1) = 1;                           % Boundary condition u_x = 0, u_0 = 1
        d(end) = 1;                         % Boundary condition u_ny = 1
        uNew = Thomas3(a,b,c,d);
        u(:,j) = uNew;                      % Updates u matrix
        uNew = uNew';
        cor1 = max(abs(uNew-uOld));
        cor2 = cor1;
    
    
    
    elseif x(j) <= 1 & x(j)>0                   % Locations along plate

        a = (1/(Re*dy^2))+v(2:ny,j)/(2*dy);    % First row [1,0...][u_1] = 0
        a(end) = 0;                             % Last row [0,... 0,1][u_ny] = 1
        b = (-2/(Re*dy^2)) - uOld/dx;
        b(1) = 1;
        b(end) = 1;
        c = (1/(Re*dy^2))-v(1:ny-1,j)/(2*dy);
        c(1) = 0;
        for k = 1:ny
            d(k) = (-uOld(k)*(u(k,j-1)))/dx;
        end
        d(1) = 0;
        d(end) = 1;
        uNew = Thomas3(a,b,c,d);
        u(:,j) = uNew;
        uNew = uNew';
        cor1 = max(abs(uNew-uOld));

        for i=2:ny                             % Solves y past plate with y_1 = 0
        vNew(i-1) = (-dy/(2*dx))*(u(i,j) - u(i,j-1)+u(i-1,j)-u(i-1,j-1)) + v(i-1,j);
        end

        v(2:ny,j) = vNew;
        cor2 = max(abs(vNew-vOld));
        
        
        
      
    else
        a = (1/(Re*dy^2))+v(2:ny,j)/(2*dy);
        a(end) = 0;
        b = (-2/(Re*dy^2)) - uOld/dx;
        b(end) = 1;
        c = (1/(Re*dy^2))-v(1:ny-1,j)/(2*dy);
        c(1) = 2/(Re*dy^2);
        for k = 1:ny
            d(k) = (-uOld(k)*(u(k,j-1)))/dx;
        end
        d(end) = 1;
        uNew = Thomas3(a,b,c,d);
        u(:,j) = uNew;
        uNew = uNew';
        cor1 = max(abs(uNew-uOld));

        for i=2:ny
        vNew(i-1) = (-dy/(2*dx))*(u(i,j) - u(i,j-1)+u(i-1,j)-u(i-1,j-1)) + v(i-1,j);
        end

        v(2:ny,j) = vNew;
        cor2 = max(abs(vNew-vOld));
        
        
        
    end     
    
      cor = max([cor1,cor2]);
end  
end

velocity = sqrt(u.^2 + v.^2);


figure(2)
subplot(3,1,q)
imagesc(x,y,velocity)
set(gca,'YDir','normal')

colorbar

xlabel('x')
ylabel('y')
title(['Velocity Field: Re = ', num2str(Re)])

x2 = [0,0.25,0.5,0.75,1.0,1+dx,2,3,4];
x2(1) = -1;
x1 = int32([0,0.25,0.5,0.75,1.0,1+dx,2,3,4]/dx)+1;
u1 = u(:,x1);

xPlate = [0:dx:1];
delta = 5 * sqrt(xPlate/Re);

for i = 1:9
    u1(:,i) = u1(:,i)+x2(i);
end
x2(1) = 0;



figure (3)              % Plots profile
subplot(3,1,q)
for i = 1:9
xline(x2(i),'--')       % Horizontal lines at each x
end
hold on

h = plot(xPlate,delta,'DisplayName', 'Boundary Thickness');

plot (u1,y)
xlabel('x (Local u)')
ylabel('u,y')
title({'Velocity and Boundary Layer Profile', ['Re = ', num2str(Re)]})
legend(h)
hold off


tau = (2/Re*dy) *(-3/2*(u(1,:)) + 2*u(2,:) - (1/2)*(u(3,:)));
plate = x>= 0 & x<=1;
tau = tau(plate);

figure(4) 
subplot(3,1,q)
plot(x(plate),tau)
xlabel('x')
ylabel('\tau')
title(['Wall Shear Stress Re = ', num2str(Re)])


end

%% Problem 2a
function Problem2a()
x = linspace(0,pi,101);
dx = x(2) - x(1);
dt = .99*(dx^2/2);                % Stability restriction
t = [0:dt:10*pi];
k = 1;
nx = length(x);
nt = length(t);

% Forward Euler
u = zeros(nx,nt);           % u(x,t)
u(1,:) = 0;                 % u(0,t) = 0
u(nx,:) = 0;                % u(pi,t) = 0
u(1:nx,1) = sin(x(1:nx));   % u(x,0) = sin x
u(1:nx,2) = u(1:nx,1);      % du/dt(0) = 0

for n = 2:nt-1                          % March in t
    for i = 2:nx-1                      % March in x 
        u(i,n+1) = ((k*dt^2/(dx^2))* (u(i-1,n) - 2*u(i,n) + u(i+1,n))) + 2*u(i,n)-u(i,n-1);
    end
end

figure(5)
subplot(2,1,1)
imagesc(t,x,u)
set(gca,'YDir','normal')
colorbar
xlabel('t')
ylabel('x')
title('Wave Equation - Forward Euler')

subplot(2,1,2)
surf(t,x,u,'LineStyle','none')
xlabel('t')
ylabel('x')
zlabel('u')
title('Wave Equation - Forward Euler')


% Backward Euler
u = zeros(nx,nt);           % u(x,t)
u(1,:) = 0;                 % u(0,t) = 0
u(nx,:) = 0;                % u(pi,t) = 0
u(1:nx,1) = sin(x(1:nx));   % u(x,0) = sin x
u(1:nx,2) = u(1:nx,1);      % du/dt(0) = 0

for n = 2:nt-1                          % March in t
a = ones(nx-1,1)*(k/dx^2);
b = (ones(nx,1)*((-2*k/dx^2))-1/dt^2);
c = a;
d = (1/dt^2) *(-2*u(:,n)+u(:,n-1));

sol = Thomas3(a,b,c,d);
u(1:nx,n+1) = sol(1:end);

end

figure(6)
subplot(2,1,1)
imagesc(t,x,u)
set(gca,'YDir','normal')
colorbar
xlabel('t')
ylabel('x')
title('Wave Equation - Backward Euler')

subplot(2,1,2)
surf(t,x,u,'LineStyle','none')
xlabel('t')
ylabel('x')
zlabel('u')
title('Wave Equation - Backward Euler')



% Crank Nicolson
u = zeros(nx,nt);           % u(x,t)
u(1,:) = 0;                 % u(0,t) = 0
u(nx,:) = 0;                % u(pi,t) = 0
u(1:nx,1) = sin(x(1:nx));   % u(x,0) = sin x
u(1:nx,2) = u(1:nx,1);      % du/dt(0) = 0

for n = 2:nt-1                          % March in t
a = ones(nx-3,1)*(k/(4*dx^2));
b = ones(nx-2,1)*(-2*k/(4*dx^2)-1/(dt^2));
c = a;

    for i = 2:nx-1                      % March in x
        d(i-1) = ((1/dt^2)*(-2*u(i,n)+u(i,n-1))) - (2*k/(4*dx^2))*(u(i-1,n)-2*u(i,n)+u(i+1,n))-(k/(4*dx^2))*(u(i-1,n-1)-2*u(i,n-1)+u(i+1,n-1));
    end
sol = Thomas3(a,b,c,d);
u(2:nx-1,n+1) = sol(1:end);

end

figure(7)
subplot(2,1,1)
imagesc(t,x,u)
set(gca,'YDir','normal')
colorbar
xlabel('t')
ylabel('x')
title('Wave Equation - Trapezoid Rule')

subplot(2,1,2)
surf(t,x,u,'LineStyle','none')
xlabel('t')
ylabel('x')
zlabel('u')
title('Wave Equation - Trapezoid Rule')



% Runge-Kutta 4
u = zeros(nx,nt);           % u(x,t)
v = zeros(nx,nt);           % v(x,t)
u(1,:) = 0;                 % u(0,t) = 0
u(nx,:) = 0;                % u(pi,t) = 0
u(1:nx,1) = sin(x(1:nx));   % u(x,0) = sin x
v(:,1) = 0;                 % du/dt(0) = 0

fu = @(uVec) [0;(uVec(nx+1:end-2)-2.*uVec(nx+2:end-1)+uVec(nx+3:end))/dx^2;0;uVec(1:nx)];

for j = 1:nt-1
    uVec = [v(:,j);u(:,j)];
    sol = RK4(fu,uVec,dt);
    v(:,j+1) = sol(1:nx);
    u(:,j+1) = sol(nx+1:end);
end




figure(8)
subplot(2,1,1)
imagesc(t,x,u)
set(gca,'YDir','normal')
colorbar
xlabel('t')
ylabel('x')
title('Wave Equation - RK4')

subplot(2,1,2)
surf(t,x,u,'LineStyle','none')
xlabel('t')
ylabel('x')
zlabel('u')
title('Wave Equation - RK4')

end

%% Problem 2b
function Problem2b()
dx = .0002;
dy = .0002;

x = [-1:dx:2];
y = [0:dy:2+dy];

nx = length(x);
ny = length(y);

M = sqrt(2);
B = sqrt(M^2 - 1);

u = ones(ny,nx);


tau = 0.1;
plate = [0-dx:dx:1+dx];                
yb = tau*sin(pi*plate);
for i = 2:length(plate)-1                   % Calculates v along plate
    v(i-1) = (yb(i+1) - yb(i-1))/(2*dx);
end
v = [0,v,0];

nPlate = length(v)-2;
before = x<0;
xBefore = x(before);
nBefore = length(xBefore);
after = x>1;
xAfter = x(after);
nAfter = length(xAfter);


% Region before plate has inlet velocity

% Calculates region starting at plate
for i=nBefore+1:nBefore+nPlate-1         
    q = i - nBefore ;                % First calculates u at the plate from v
    u(1,i+1) = ((dx^2)/(-B^2)) * ((v(q+2)-v(q))/(2*dx*dy) - (u(2,i)-u(1,i))/(dy^2)) + 2*(u(1,i))-u(1,i-1);
    for j = 2:ny-1                      % Solves u from B.C. at plate
        u(j,i+1) = dx^2 * (((2*u(j,i) - u(j,i-1))/dx^2)+(u(j-1,i)-2*u(j,i)+u(j+1,i))/(B^2*dy^2));
    end
    
end

for i = nPlate:nx-1
    for j = 2:ny-1                      % Solves u after plate
        u(j,i+1) = dx^2 * (((2*u(j,i) - u(j,i-1))/dx^2)+(u(j-1,i)-2*u(j,i)+u(j+1,i))/(B^2*dy^2));
    end
end

figure(9)
xlim([-1,2])
ylim([-2,2])
hold on
imagesc(x,y,u)

imagesc(x,-y,u)
set(gca,'YDir','normal')
plot(plate,yb,'r')
plot(plate,-yb,'r')
colorbar
xlabel('x')
ylabel('y')
title('Streamwise Velocity')
hold off

figure(10)
cp = 2*(1-u(1,nBefore+1:nBefore+nPlate));
p = plate(2:end-1);
plot(p,cp)
ylabel('C_p')
xlabel('Chord x/c')
title('Pressure Distribution')

end
%% Tridiagonal Thomas Algorithm Function
function xBar = Thomas3(a,b,c,d)
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

%% Runge Kutta 4 
function y = RK4(f,v,dt)
    y = v;

        k1 = dt*f(v);
        k2 = dt*f(v+k1/2);
        k3 = dt*f(v+k2/2);
        k4 = dt*f(v+k3);
        y = v + (1/6)*(k1+2*k2+2*k3+k4);

end

