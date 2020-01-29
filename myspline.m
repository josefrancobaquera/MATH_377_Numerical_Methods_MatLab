% HW7 - Jose Franco Baquera 
% MATH 377
% April 2, 2019

function pp = myspline(x,y)

% This function finds the natural cubic spline given data points.

% Assumptions:
% x are nodes (x value) = [x0,x1,...,xn]
% y are values (y value) = [y0,y1,...,yn]

% Output pp is a structure.

% Note: pp.coef is (length(x)-1) by 4 matrix. One row stores the 4
% coefficients of a cubic piece; degree 3 first.
% Coefficients are local: if coef = [a,b,c,d], then first piece 
% is a(x-x0)^3+b(x-x0)^2+c(x-x0)+d

% Check if input is correct first.

% Check that both inputs are vectors. A vector is an array that has a 
% size of 1-by-N or N-by-1, where N is a nonnegative integer. In other
% words, it checks that inputs are not a matrix.
if ~isvector(x)
    error('x needs to be a vector')
end
if ~isvector(y)
    error('y needs to be a vector')
end

% Calculate n. There are n+1 points.
n = length(x) - 1;  

% Make sure that there are at least 3 nodes and that x and y have
% the same length.
if n<2
    error('There needs to be at least 3 nodes.')
end
if length(y)~= n+1
    error('x and y should have the same length.')
end

% Turn both x and y into column vectors.
x = reshape(x,[n+1,1]); 
y = reshape(y,[n+1,1]);

% Calculate the h's. It is important to note that h's 
% start at index 1 in both MATLAB and in the equation and 
% go all the way to n. Note: the abs is not needed 
% but is included just in case.
h =  abs(x(2:n+1) - x(1:n));

% Step 1 - Construct the tridiagonal matrix:

% It is important to note that the tridiagonal 
% matrix will be a n-1 by n-1 square matrix.

% First, compute the diagonal entries of the matrix.
diagonalEntries = (h(1:n-1) + h(2:n))/3;

% Compute the actual tridiagonal matrix. The function 'gallery' returns 
% the coordinates of all the non-zero entries of the corresponding matrix.
% The 'full' function contructs the actual matrix with both nonzero and 
% zero entries. Assign result to variable. NOTE: If there are 3 nodes
% only, treat the matrix a special way (i.e. matrix will be 1 by 1). Else, 
% continue normally. 
if n == 2
   tridiagonalMatrix = diagonalEntries;
else
    tridiagonalMatrix = full(gallery('tridiag', n-1, h(2:n-1)/6, diagonalEntries, h(2:n-1)/6));
end
    
% Step 2 - Construct the right hand side:

% We construct the right hand side by doing vector subtraction/division and using
% the formula provided to us on the notes.
rightHandSide = ((y(3:n+1)-y(2:n))./h(2:n)) - ((y(2:n)-y(1:n-1))./h(1:n-1));

% Step 3 - Solve for z1,...,zn-1:

% We note that there are a total of n+1 z's.
% Since this is a natural cubic spline, z0=zn=0.
% However, since indeces start at 1 in MATLAB, 
% z0=z1 and zn=zn+1. First, we make a vector of the 
% correct size with all zeros. 
z = zeros(n+1, 1);

% Now we solve for z1,...,zn-1. However, since indeces start 
% at 1 in MATLAB, we will be storing z1,...,zn-1 in
% z2,...,zn.

z(2:n) = tridiagonalMatrix\rightHandSide;

% Step 4 - Compute coefficients:

% We will have n number of cubic pieces. That is, one row stores the 4
% coefficients of a cubic piece; degree 3 first.
% First create a generic matrix will all zeros. The final coefficients will
% be stored here.
coefs = zeros(n,4);

% Note: ith row is the coefs for qi.

% Compute the coefficients for degree 3. 
coefs(:,1) = (z(2:n+1) - z(1:n))./(6*h(1:n));

% Compute the coefficients for degree 2. 
coefs(:,2) = z(1:n)*(.5);

% Compute the coefficients for degree 1.
coefs(:,3) =  ((y(2:n+1) - y(1:n))./h) - (h.*(z(2:n+1)/6 +z(1:n)/3));

% Compute the coefficients for degree 0.
coefs(:,4) =  y(1:n);

% pp is a structure. Assign the corresponding instance variables to their 
% appropriate values.
pp.form = 'pp';
pp.breaks = reshape(x,[1,n+1]);
pp.coefs = coefs;
pp.pieces = n;
pp.order = 4;
pp.dim = 1;


    
