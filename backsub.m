
%% HW5 - Jose Franco Baquera 

% February 28th, 2019
% File Name: backsub.m
% MATH 377

function x = backsub(U, b)

% Write a function x = backsub(U, b) to do backward substitution from scratch. The input
% U must be an invertible upper triangular matrix (of any size) and b is a column vector of
% appropriate length. The function returns the solution of Ux=b.

[m, n] = size(U); % Retrive the size (row by column) of the matrix.
[l, p] = size(b); % Retrive the size of the vector. 

% First check if input is correct. 

% Check if matrix U is n by n. (NOTE: Invertible matrices can only be square.)
if m~=n 
    error('The matrix is not invertible (i.e. the first input needs to be a square).')
end

% Check if vector b is n by 1.
if l~=n || p~=1 
    error('The second input does not have appropriate size (i.e. must be n by 1).')
end

% Check if matrix U is upper triangular. 
if istriu(U) == 0
   error('The matrix is not upper triangular. Try again.')    
end

% Check if the matrix is invertible by checking if there is a pivot in every column.  
for i = 1:n % Go through each column.
    
   if U(i, i) == 0 % If no pivot is found, matrix is not invertible. Throw an error.
      error('The matrix is not invertible (i.e. does not contain n pivots).')
   end % end if
   
end % end for

% If we get to this place, we assume that the input is all correct. 
% Now we need to find the solution to Ux=b using one for loop.
x = ones(n, 1); % Create a vector the same size of b with all one's.

% For loop that will perform "backward substitution". We need to go from
% one row to another, starting with the last row.
for i = m:-1:1 % For loop that will count "backwards".
    
    if i == m % Treat the last row differently than the rest. 
       x(m) = b(m)/U(m,m); % Solve for the "last" variable and store it in x.
    else
       tempVar = b(i) - dot( U(i, i+1:n)', x(i+1:m) ); % Do dot product to find the "top" part.
                                                       % of the division templete provided
                                                       % to us in class.
       x(i) = tempVar /  U(i,i); % Solve for the "ith" variable and store it in x.  
    end % end if else
    
end % end for.






    