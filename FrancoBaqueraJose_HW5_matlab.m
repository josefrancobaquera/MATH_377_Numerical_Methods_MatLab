%% HW5 - Jose Franco Baquera 

% February 28th, 2019
% File Name: FrancoBaqueraJose_HW5_matlab.m
% MATH 377

% format short to display matrices better.
format short;
%% Ch3-22: Basic Matrix Operations

% (a) Create a random 4 × 4 matrix A (each entry follows standard normal distribution)
disp('Random Matrix A:'); % Display a meaningful message.
A = randn(4,4) % This creates a 4 by 4 matrix called A with random entries from the standard normal 
               % distribution. 
               
% (b) Display the 3rd column of A.
disp('Column 3 of Matrix A:'); % Display a meaningful message.
A(:,3) % This displays the 3rd column of A.

% (c) Display the upper right 2 × 2 submatrix of A.
disp('Upper Right 2 x 2 Submatrix of A:'); % Display a meaningful message.
A(1:2,3:4) % This displays the upper right 2 x 2 submatrix.  

% (d) Compute the sum of all entries of A.
disp('Sum of All Entries of A:'); % Display a meaningful message.
sum(A(:)) % Computes the sum of all rows and columns of matrix A. 

% (e) Define B to be the transpose of A.
disp('Let B Equal the Transpose of A. B is Equal to:'); % Display a meaningful message.
B = A' % The ' character transposes a matrix. 

%% Ch3-23 (How powerful is your computer (or the lab computer)?)

% format long to make calculations more precise. 
format long;

% (a) Create a 5000 by 5000 random matrix R, where each entry of R is standard normal
% distribution. Define r to be the product of R and a 5000 × 1 vector of all 1's.
R = randn(5000,5000); % This creates a random 5000 by 5000 matrix called R with standard entries.
allOnesVector = ones(5000,1); % Creates a 5000 by 1 vector is all one's.
r = R*allOnesVector; % Define r to be the product of R and allOnesVector.

% (b) Measure how much time does it take to solve Rx=r using Matlab's built-in function. You
% may want to run several times so that the time stabilizes. Use denseGE to store the
% elapsed time.
tic; % Start timer.
x = R\r; % Solve Rx=r using Matlab's built-in function. Store it in x and supress printing.
denseGE = toc; % Stop timer and store elapsed time to variable denseGE.
disp('Total Elapsed Time (Seconds):'); % Display a meaningful message.
denseGE % Display the elapsed time.

% (c) Calculate how many gigaflops your computer can do in a second. (This is related to
% Problem 17.) Display your answer using the command fprintf.
% The stabilized time was around 1.3288 seconds. 
% We note that to solve a 5000 by 5000 system, it takes roughly about (2/3)n^3 flops, or about
% (2/3)(5000)^3 flops.
totalNumberOfFlops = (2/3)*5000^3;
% However, since my computer computed this amount of flops in about 1.3288
% seconds, we still have to find the amount of flops it executed in one
% second. We use denseGE to compute the most accurate number of flops
% my computer can do. That is:
numberOfFlopsPerSecond = totalNumberOfFlops/denseGE;
% NOTE: One gigaflop equals to 10^9 flops. 
fprintf( 'Total Amount of Flops Per Second (in GIGAFLOPS): %.2f', numberOfFlopsPerSecond/(10^9) )

%% Ch3-24

% There is usually no good reason for ever computing the inverse of a matrix in practice. Take
% the R and r from Problem 23, but solve Rx=r by first finding the inverse of R, and then
% performing the multiplication R^-1r. Measure how much time it takes. Use denseINV to store
% the elapsed time. (According to the top of page 20, you can expect denseINV to be at least
% twice as much as denseGE.)
tic; % Starting timer.
inverseR = inv(R); % Compute the inverse.
x = inverseR*r; % Solve the system and store result.
denseINV = toc; % Stop timer and store elapsed time to variable denseINV.
disp('Total Elapsed Time (Seconds):'); % Display a meaningful message.
denseINV % Display the elapsed time.
% We note that denseINV is at least twice as much as denseGE.

%% Ch3-25

% Write a function x = backsub(U, b) to do backward substitution from scratch. The input
% U must be an invertible upper triangular matrix (of any size) and b is a column vector of
% appropriate length. The function returns the solution of Ux=b.

% Please look at backsub.m file.

%% Ch3-26 (More time comparison)

% (a) Define U to be the upper triangular portion of a 5000 × 5000 matrix full of 1's. Define b
% to be the product of U and a 5000 × 1 vector of all 1's.
U = ones(5000, 5000); % Defines a 5000 by 5000 matrix with all ones.
U = triu(U); % Make U upper triangular. 
allOnesVector = ones(5000,1); % Creates a 5000 by 1 vector is all one's.
b = U*allOnesVector; % Defines b to be the product of U and a 5000 by 1 vector of all 1's. 

% (b) Solve Ux=b using Matlab's built-in function. Use upperMatlab to store the elapsed time.
tic; % Starting timer.
x = U\b; % Solve the system and store the result.
upperMatlab = toc; % Stop timer and store elapsed time to variable upperMatlab.

% (c) Solve Ux=b using your own function backsub. Use upperMine to store the elapsed time.
% Check that the solution is the same as the one from (b) to make sure your code is correct.
tic; % Starting timer.
backsubResult = backsub(U, b); % Solve the system and store the result.
upperMine = toc; % Stop timer and store elapsed time to variable upperMine.

% Check if both results are equal.
if x == backsubResult
   disp('backsub function works! Yay! :)'); % Display a meaningful message.
else
   disp('backsub function does not work! :('); % Display a meaningful message.
end

% (d) Finally display all 4 times (denseGE, denseINV, upperMatlab, upperMine) using table.
% Comment on them.
TypeOfTime = { 'denseGE'; 'denseINV'; 'upperMatlab'; 'upperMine' }; % Create the first column.
TimeInSeconds = [ denseGE, denseINV, upperMatlab, upperMine ]'; % Create the second column.
myTable = table ( TypeOfTime, TimeInSeconds ) % Display all the 4 times using a table.

% Comments on Table: We note that matlab's built-in function to solve 
% linear systems was a little bit more than twice as fast as first finding 
% the inverse and multiplying it with r. To solve Ax=b without first finding 
% the inverse takes 2/3n^3 flops. However, if we wanted to take the inverse first it 
% would take 2n^3 flops, which is about 3 times more flops then not finding the
% inverse. Because of this, finding the solution by using the matrix's
% inverse takes more time than using matlab's built-in function. The time
% it took matlab to solve an upper triangular system using its built-in
% function was much faster than our backsub function implementation. This is 
% probably because our function performs more flops than the built-in matlab function.
% We conclude that the matlab built-in function is more effecient than our implementation
% when solving upper triangualar systems. Futhermore, the time it takes to solve an upper
% triangular matrix is much more faster than a dense one. Matlab's built-in
% function probably first checks if the matrix is upper triangular first, then continues 
% to solve the system using only backwards substitution, which only takes O(n^2) flops.

