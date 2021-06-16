%Chris Hopp
%915866326
%ENG-180 Project 1: Gauss Elimination
%10/6/2020


clc
clear all
%#ok<*SAGROW>
%#ok<*NOPTS>

%% Problem 1 Input
for j=1:3
    n = 10^j;                   % Iterates through n=10,100,1000
    a = -1*ones(n-1,1);         % Subdiagonal
    b = 2*ones(n,1);            % Diagonal
    c = a;                      % Superdiagonal
    d = zeros(n,1);             % RHS of equation
    d(1) = 1;
    d(n) = 1;

    xVal = THOMAS3(a,b,c,d);    % x-vector solution
    x{j} = xVal;                % Data for solution table
    maxVal(j) = max(xVal);
    minVal(j) = min(xVal);
    nValue(j) = n;
end

Problem1_X = cell2mat(x(1));     % Problem 1 solution
Problem1 = table(nValue', maxVal', minVal', x', 'VariableNames', {'n','Max', 'Min', 'x'});


%% Problem 2 Input
for j=1:3
    n = 10^j;                   % Iterates through n=10,100,1000
    N = n/2;
    a = cell(N-1,1);            % Allocation for diagonals and RHS
    b = cell(N,1);             
    c = cell(N-1,1);            
    d = cell(N,1);              
    
    for i=1:N-1                 % Forms sub and superdiagonals
        a{i} = eye(2);
        c{i} = eye(2);
    end
    
    for i=1:N                   % Forms diagonal and RHS
        b{i} = [-2,-1;0,-2];
        d{i} = [-1;0];
    end
    d{1} = [-2;-1];
    d{N} = [-2;-1];
    
    xVal = THOMAS3_BLOCK(a,b,c,d);    % x-vector solution
    x{j} = xVal;                % Data for solution table
    maxVal(j) = max(xVal);
    minVal(j) = min(xVal);
    nValue(j) = n;
end

Problem2_X = cell2mat(x(1));     % Problem 2 solution
Problem2 = table(nValue', maxVal', minVal', x', 'VariableNames', {'n','Max', 'Min', 'x'});


%% Problem 3a Input
for j=1:3                           
    n = 10^j;                    % Iterates through n=10,100,1000
    a = ones(n-2,1);             % Lower subdiagonal
    b = -4*ones(n-1,1);          % Subdiagonal
    c = 6*ones(n,1);             % Diagonal
    d = b;                       % Superdiagonal
    e = a;                       % Upper superdiagonal

    f = zeros(n,1);              % RHS of equation
    f(1)=3; f(2)=-1; f(n-1)=-1; f(n)=3;

    xVal = THOMAS5(a,b,c,d,e,f);    % x-vector solution
    x{j} = xVal;                % Data for solution table
    maxVal(j) = max(xVal);
    minVal(j) = min(xVal);
    nValue(j) = n;
end

Problem3a_X = cell2mat(x(1));     % Problem 3a solution
Problem3a = table(nValue', maxVal', minVal', x', 'VariableNames', {'n','Max', 'Min', 'x'});


%% Problem 3b Input
for j=1:3
    n = 10^j;                       % Iterates through n=10,100,1000
    N = n/2;                        % Length of cell array holding 2x2 blocks
    a = cell(N-1,1);                % Allocation for diagonals and RHS
    b = cell(N,1);
    c = cell(N-1,1);
    d = cell(N,1);


    for i=1:N                       % Forms diagonal
        b{i} = [6,-4;-4,6];
    end
    
    for i=1:N-1                     % Forms sub and superdiagonal
        a{i} = [1,-4;0,1];
        c{i} = [1,0;-4,1];
    end
    
    for i=2:N-1                     % Forms RHS
        d{i} = [0;0];
    end
    d{1} = [3;-1];
    d{N} = [-1;3];

    xVal = THOMAS3_BLOCK(a,b,c,d);    % x-vector solution
    x{j} = xVal;                % Data for solution table
    maxVal(j) = max(xVal);
    minVal(j) = min(xVal);
    nValue(j) = n;

end

Problem3b_X = cell2mat(x(1));     % Problem 3b solution
Problem3b = table(nValue', maxVal', minVal', x', 'VariableNames', {'n','Max', 'Min', 'x'});

%% Problem 4 Input
for j=1:3
    n = 10^j;                   % Iterates through n=10,100,1000
    a = -1*ones(n-1,1);         % Subdiagonal
    b = 2*ones(n,1);            % Diagonal
    c = a;                      % Superdiagonal
    d = zeros(n,1);             % RHS of equation
    d(1) = 1;
    d(n) = 1;

    xVal = LU(a,b,c,d);    % x-vector solution
    x{j} = xVal;                % Data for solution table
    maxVal(j) = max(xVal);
    minVal(j) = min(xVal);
    nValue(j) = n;
end

Problem4_X = cell2mat(x(1));     % Problem 1 solution
Problem4 = table(nValue', maxVal', minVal', x', 'VariableNames', {'n','Max', 'Min', 'x'});


%% Solution Output
Problem1_X 
Problem1
Problem2_X
Problem2
Problem3a_X
Problem3a
Problem3b_X
Problem3b
Problem4_X
Problem4



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

%% Tridiagonal LU Factorization Function
function xBar = LU(a,b,c,d)

    n = length(b);                      % Length of total diagonal
    u = zeros(n,1);                     % Vector for U matrix diagonal. U Matrix superdiagonal = c.  All other elements 0.
    l = zeros(n,1);                     % Vector for L matrix subdiagonal.  L diagonal is ones.  All other elements 0.
    a = [0;a];                          % Pads a vector
    y = zeros(n,1);                     % Allocates for y vector Ly=d
    xBar = zeros(n,1);                  % Allocates for x vector UxBar=y
    
    u(1) = b(1);                        % Forms u and l vectors to produce A = LU
    for i=2:n
        l(i)=a(i)/u(i-1);
        u(i)=b(i)-l(i)*c(i-1);
    end
    
    y(1) = d(1);                        % Solves for y vector by y=L/d solving downwards
    for i=2:n
        y(i) = d(i)-(l(i)*y(i-1));
    end
    
    xBar(n)=y(n)/u(n);                  % Solves for x vector Ux=y solving upwards
    for i=n-1:-1:1
        xBar(i)=(y(i)-c(i)*xBar(i+1))/u(i);
    end
    xBar = xBar';
end

%% 2x2 Inverse Function
function invMat = INVERSE(A)
        detA = A(1,1)*A(2,2) - A(1,2)*A(2,1);               % Calculates 2x2 determinant
        invMat = (1/detA)*[A(2,2),-A(1,2);-A(2,1),A(1,1)];  % Calculates inverse by 1/det[d,-b;-c,a]
end