%% HW10 - Jose Franco Baquera 

% February 2nd, 2019
% File Name: FrancoBaqueraJose_HW10_Matlab.pdf
% MATH 377

% Format short in order to have easier outputs.
format short;

%% Chapter 6: Problem 20

% First create a random 10×5 Gaussian matrix A and a random 10×1 vector b. It is likely that
% Ax = b has no solution. Use Matlab's built-in command to find its least squares solution,
% name it x_matlab. Now use Cholesky factorization (see the table on page 50) to find its least
% squares solution, name it x_chol (you can use Matlab's function chol). Finally check that
% these two solutions are the same/close by computing their distance.

% Create a random 10x5 Gaussian matrix A.
A = randn(10, 5);
% Create a random 10x1 vector b.
b = randn(10, 1);

% Least-squares solution using Matlab's built-in function.
x_matlab = A\b;

% Use Cholesky factorization to find solution.

% Step 1: Form B = A'A, f = A'b
B = A'*A;
f = A'*b;

% Step 2: Find Cholesky factorization B = LL'
LT = chol(B);
L = LT';

% Step 3: Solve Ly = f
y = L\f;

% Step 4: Solve L'x = y
x_chol = LT\y;

% Checking that check these two solutions are the same/close by computing their distance.
distanceVer1 = norm(x_matlab);
distanceVer2 = norm (x_chol);

% Checking the difference of these two solutions. The difference should be
% "small".
differenceBetweenSolutions = x_matlab - x_chol;

% Print a meaningful message.
fprintf ( "The norm of x_matlab is: %f\n", distanceVer1 );
fprintf ( "The norm of x_chol is: %f\n", distanceVer2 );
fprintf ( "The difference in values for each entry is the following: \n" );
differenceBetweenSolutions
fprintf ( "We conclude that both methods return extremly similar answers. \n\n" );

%% Chapter 6: Problem 21

% Use the data of the following 14 schools to predict other college's scores.

% (a) First use the model y = w0 + w1x1 + w2x2 + ..... + w8x8 to fit the
% above data and find weights. 

fprintf ( "PART A\n");

% First, create matrix A that has all x's values.
A = [1 17 71 4 7 93 1 97 58;...
     1 13 70 2 9 94 1 98 57;...
     1 15 74 2 8 93 4 97 46;...
     1 18 68 1 9 94 11 96 55;...
     1 14 70 1 8 94 1 98 43;...
     1 16 68 1 10 93 6 96 50;...
     1 31 69 1 8 93 14 95 46;...
     1 31 65 1 9 97 4 97 58;...
     1 25 79 1 8 94 6 96 44;...
     1 14 86 2 9 94 11 96 43;...
     1 23 68 0.3 8 95 6 97 33;...
     1 28 69 0 11 99 6 96 53;...
     1 22 67 2 8 97 21 98 33;...
     1 7 61 0 9 94 25 97 21];

% Now create the y vector (i.e. the corresponding scores).
y = [100 98 96 94 94 93 93 92 91 90 90 89 89 88]';

% Firstly, compute "left-hand" side.
LHAND = A'*A;

% Secondly, compute "right-hand" side.
RHAND = A'*y;

% Now we can find all 9 weights. 
WEIGHTS = LHAND\RHAND;

% Print a meaningful message.
fprintf ( "The 9 weights are the following:\n " );
WEIGHTS

fprintf ( "The equivalent funtion is the following:\n " );
fprintf ( "y = " );
fprintf ( "%f + (%f)x1 + (%f)x2 + (%f)x3 + (%f)x4 + (%f)x5 + (%f)x6 + (%f)x7 + (%f)x8\n\n",...
        WEIGHTS(1), WEIGHTS(2), WEIGHTS(3), WEIGHTS(4), WEIGHTS(5), WEIGHTS(6), ...
        WEIGHTS(7), WEIGHTS(8), WEIGHTS(9));

% (b) Second estimate the scores of the 4 colleges below.

fprintf ( "PART B\n");

% Create matrix AESTIMATE of all x's values for the colleges we want to
% estimate. 
AESTIMATE = [ 1 18 74 0.2 9 91 14 94 46;...
              1 27 74 1 9 94 21 95 47;...
              1 24 68 5 9 97 6 96 49;...
              1 29 69 2 10 93 14 95 41];
  
% Estimate the scores of the 4 collages below:
yESTIMATE = AESTIMATE*WEIGHTS;

% Print a meaningful message.
fprintf ( "The score estimate for Wash. & Lee college is %f.\n", yESTIMATE(1));
fprintf ( "The score estimate for Hamilton college is %f.\n", yESTIMATE(2));
fprintf ( "The score estimate for Welsleyan college is %f.\n", yESTIMATE(3));
fprintf ( "The score estimate for Colby college is %f.\n", yESTIMATE(4));

%% Chapter 6: Problem 22

% Line fitting using both least squares and SVD on (0,1), (1,3), (2,4), (3,4).

% (a) Find the equation of the best fitted line y = ax + b using least squares.

% y = ax + b => y = b + ax
% First, create the required matrices.
LSMA = [1 0;...
       1 1;...
       1 2;...
       1 3];
LSMB = [1 3 4 4]'

% Second, find the coefficients.
PDLSMA = LSMA'*LSMA;
TEMPLSMP = LSMA'*LSMB;
LSC = PDLSMA\TEMPLSMP;

% Print a meaningful message.
fprintf ( "PART A\n " );
fprintf ( "a = %f and b = %f\n ", LSC(2), LSC(1) ); 

fprintf ( "Therefore, the equivalent funtion is the following:\n " );
fprintf ( "y = %f + (%f)x\n\n",  LSC(2), LSC(1));

% (b) Find the equation of the best fitted line y = cx + d using SVD.

% First find D (i.e. each point is a row of D).
D =  [0 1; 1 3; 2 4; 3 4];

% Second find E.
C = mean(D);
E = D - C;
[E1, E2, E3] = svd(E);
% Compute the slope of line.
slope = E3(2)/E3(1);

% Therefore, the line is equal to y = slope(x - c1) + c2
% Print a meaningful message.
fprintf ( "PART B\n " );

fprintf ( "The best fitted line is the following:\n " );
fprintf ( "y = %f(x-%f) + %f\n", slope, C(1), C(2));
fprintf ( "This equation can be written (equally) as the following:\n " );
fprintf ( "y = %fx + %f\n\n", slope, slope*-1*C(1)+C(2));

% (c) Plot these 4 points and both lines in the same graph. Make sure you have a legend.

% Define both equations to graph. 
yOne = @(x) 1 + 1.5*x;
yTwo = @(x) slope*x + slope*-1*C(1)+C(2);

% Plot both line fittings, as well as the four points.

fprintf ( "PART C\n " );
figure(1)
xIntervalStep = 0:0.01:3;
xPoints = [0 1 2 3];
yPoints = [1 3 4 4]; 
plot(xIntervalStep, yOne(xIntervalStep), xIntervalStep, yTwo(xIntervalStep), xPoints, yPoints, 'ko') %, xIntervalStep, yTwo);
legend('Best Line Using Least Squares','Best Line Using SVD');
title('Chapter 6: Problem 22 - Part C')
