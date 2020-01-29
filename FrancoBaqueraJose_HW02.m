%% HW2 - Jose Franco Baquera (Introduction)

% February 7th, 2018
% File Name: FrancoBaqueraJose_HW02.m
% MATH 377

%% Ch1-11

% (a) Define a row vector v whose entries are [1, 0.01, 0.02, ...  , 5].  Do not print.
v = 1:0.01:5; % This defines a row vector with entries [1, 1.01, 1.02, ... , 5]. 
              % It includes a semicolon so that it is not printed. 
                   
% (b) Define a 10×10 matrix whose entries are the integers 1 to 100 ordered column-wise. Do not print.
matrix = [1:10; 11:20; 21:30; 31:40; 41:50; 51:60; 61:70; 71:80; 81:90; 91:100 ]';
   % This defines a 10x10 matrix whose entries are integers 1 to 100 
   % ordered column-wise. It includes semicolon so that it is not printed. 
   % The ' character is included since it "transposes" the matrix. 

% (c) Plot y = cosx and y = x on the interval [0,1] in the same graph. Submit plot as .png.
x = 0:0.01:1; % Variable that will store the interval [0,1]. We use 0.01 increments so that 
              % the graph looks "precise." Do not print it this particular variable.
              
figure(1) % This graph will be labeled as figure 1.

plot(x, cos(x), 'k', x, x, 'r') % Plot y = cos x and y = x on the same graph and on interval [0,1].
                                % y = cos x will be black in color and y = x will be red in color.     
                                
title( 'y = cos(x) (black) and y = x (red) on interval [0,1]' ) % Label the title of the graph since
                                                                % it is good practice to do so. 

%% Ch1-12

% Given the equation x^2 - 100000x + 1 = 0, use formula (1.5) and formula (1.6) to compute
% two different approximations of x1 in Matlab. Print both values in the long format. Which
% formula produces a better approximation of x1?

format long; % This prints numbers in long format. 

disp('Value of x1 using equation (1.5)') % Display a meaningful message.
x1a = (100000 + sqrt((-100000)^2 - 4*1*1))/2; % Value of x1 using equation 1.5. 
x1a % Print the value of x1 that was computed by formula 1.5. 

disp('Value of x1 using equation (1.6)') % Display a meaningful message.
x1b = (-2*1)/(-100000 + sqrt((-100000)^2 - 4*1*1)); % Value of x1 using equation 1.6. 
x1b % Print the value of x1 that was computed by formula 1.6. 

% Conclusion: We note that formula 1.6 produced a better approximation of x1.
disp('Conlusion: Formula (1.6) produced a better approximation of x1.') % Display a meaningful message.

%% Ch1-13

% Given p(x) = x^9 -18x^8 + 144x^7 - 672x^6 + 2016x^5 - 4032x^4 + 5376x^3 - 4608x^2 + 2304x - 512  

% (a) Let n(x) be the nested form of p(x). Define p(x), n(x), and s(x) = (x - 2)^9 in Matlab.
% All 3 functions need to be take in vectors.

% NOTE: The dot in certain places allows all the functions to take in vecotrs.  
% This dot can ONLY be used in the following way: '.*', './', and '.^'
p = @(x) x.^9 - 18*x.^8 + 144*x.^7 - 672*x.^6 + 2016*x.^5 - 4032*x.^4 + 5376*x.^3 - 4608*x.^2 + 2304*x - 512;
   % p is the original equation. The extra dots after certain x's allows the function to take in
   % vectors as input. It is important to note that the last x does not require a dot. 
n = @(x) ((((((((x-18).*x + 144).*x - 672).*x + 2016).*x -4032).*x + 5376).*x -4608 ).*x + 2304).*x -512;
   % n is the nested equation. The dots before all the multiplication signs allow the equation to take in
   % vectors as input. 
s = @(x) (x-2).^9;
   % This function was given on the assignment. The extra dot after the left parathensis allows the 
   % function to take in vectors as input. 

% (b) p(x) = n(x) = s(x) algebraically, but makes a difference in computers. Plot all three
% functions on the interval [1.92, 2.08] (with 0.001 increment). Use the command subplot
% to place all 3 plots in a 1 × 3 grid. Be sure to give a title to each subplot to indicate
% which is which.

x = 1.92:0.001:2.08; % Variable that will store the interval [1.92,2.08] with 0.001 increment.
figure(2) % Figure 2 will be these three subplots.

subplot(1,3,1) % Assign the first subplot to equation p(x).
plot(x,p(x)) % Plot p(x).
title({'Graph of Equation p(x)',' '}) % Give a meaningful title.

subplot(1,3,2) % Assign the second subplot to equation n(x).
plot(x,n(x)) % Plot n(x).
title({'Graph of Equation n(x)', ' '}) % Give a meaningful title.

subplot(1,3,3) % Assign the third subplot to equation s(x).
plot(x,s(x)) % Plot s(x).
title({'Graph of Equation s(x)', ' '}) % Give a meaningful title.

% (c) Explain the difference in three plots in (b). 
% The difference in the three plots in (b) is quite drastic. The graph for equation p(x)
% has a lot of oscillation as the graph approaches the root. In other words, there is are lot
% of oscillation as we approach p(x) = 0. The graph for equation n(x) also has this 
% oscillation as we approach n(x) = 0, but it has "less oscillation" than p(x). The graph 
% of s(x) has no oscillation as we approach s(x) = 0. To summarize, the amount of 
% oscillation as we approach the root decrease from p(x) to s(x). The reason why 
% this is the case is because p(x) has the greatest number of FLOPS while s(x)
% has the least number of FLOPS. n(x) has less FLOPS than p(x), but more FLOPS than
% s(x). These three graphs demonstrate that s(x) has a better approximation of the root
% while p(x) has the worst. We conclude that more oscillation correlated directly to 
% less precise approximation, which also directly correlates to more FLOPS.

disp( 'The difference in the three plots in (b) is quite drastic. The graph for equation p(x)' )
disp( 'has a lot of oscillation as the graph approaches the root. In other words, there is are lot')
disp( 'of oscillation as we approach p(x) = 0. The graph for equation n(x) also has this' )
disp( 'oscillation as we approach n(x) = 0, but it has "less oscillation" than p(x). The graph' ) 
disp( 'of s(x) has no oscillation as we approach s(x) = 0. To summarize, the amount of' )
disp( 'oscillation as we approach the root decrease from p(x) to s(x). The reason why' )
disp( 'this is the case is because p(x) has the greatest number of FLOPS while s(x)' )
disp( 'has the least number of FLOPS. n(x) has less FLOPS than p(x), but more FLOPS than' )
disp( 's(x). These three graphs demonstrate that s(x) has a better approximation of the root' )
disp( 'while p(x) has the worst. We conclude that more oscillation correlated directly to' ) 
disp ( 'less precise approximation, which also directly correlates to more FLOPS.' )

%% Ch2-1

% Bisection and Newton comparison.

% (a) Use the Bisection method to do Example 2.1, with the stopping criteria |pn+1 - pn| <= 1e-5, 
% and initial interval a0 = 0, a1 = 1. How many iterations was run? (meaning what is
% the n value when |pn+1 - pn| <= 1e-5 is reached?) How is this compared to the theoretical
% value in (2.1)?

clear x; % Clear x from previous problem. 
syms x; % Make x a symbolic variable. 
f = @(x) (0.5)^x - x; % Define the original equation used in example 2.1.

N = 20; % Set a max number of iterations in case we never reach the break case.  
p = ones(N,1); % Make a vector with N number of ones. The answers 
                             % of every iteration of the bisection method
                             % will be stored in this vector.
                             
tol = 1e-5; % Same as 10^(-5). This is the tolerance level provided to us by the problem. 

% Set the initial interval values.
% We note that f(0) > 0 and f(1) < 0.
a = 0;
b = 1;

pn = (a+b)/2; % Computing pn. NOTE: p0 does not count as an iteration.

numOfIterations = 0; % Define a new variable to keep track of the number of iterations.

% For loop that will go until N iterations or until the breaking case is reached.
for i = 1:N
   
   if f(pn) < 0  % If p0 < 0, update "b" since "b" will always point to negative values. 
      b = pn;    
   else
      a = pn; % If p0 > 0, update "a" since "a" will always point to positive values.  
   end
   
   if abs( (a+b)/2 - pn ) <= tol % Check the stoping criteria/tolerance level.
      break; % Break if the tolerence level was met.     
   end
   
   % Continue with the iteration process if the tolerance level was not met. 
   p(i) = (a+b)/2; % Split interval and choose middle vlaue. Store it in the array.
   pn = (a+b)/2; % Update pn.
   numOfIterations = numOfIterations + 1; % Increment the total number of iterations.
   
end

disp( 'Here are the 15 iterations for the Bisection Method. NOTE: P0 is NOT shown since it is not part of the iteration process.' ) 
% Printing the list of numbers that got closer to the real root. 
p(1:15)

% How many iterations was run? (meaning what was the n value?).
% The total number of iterations was 15.
disp( 'How many iterations was run? (meaning what was the n value?)' ) 
disp( '(Note: We do not count p0 as an iteration)' )
disp( 'The total number of iterations was the following:' ) 

numOfIterations 

% How is this compared to the theoretical value in (2.1)?
% The theoretical value computed in 2.1 was 15 (rounded).
% Therefore, n should be greater than or equal to 15. For our particular problem,
% we got the total number of iterations equal to 15. Because 15 >= 15, then
% the equation does check out. 

disp( 'How is this compared to the theoretical value in (2.1)?' ) 
disp( 'The theoretical value computed in 2.1 was 15 (rounded).' ) 
disp( 'Therefore, n should be greater than or equal to 15. For our particular problem,' )
disp( 'we got the total number of iterations equal to 15. Because 15 >= 15, then' )
disp( 'the equation does check out.' )

% (b) Use the Newton's method to do Example 2.1, with the stopping criteria |xn+1 - xn| <= 1e-5, and initial value x0 = 1. 
% How many iterations was run?

% We note, that we can also use function f that was previously declared.
fderi = @(x) -(1/2)^x*log(2) - 1; % Declare the FIRST derivate of function f.

N = 20; % Set a max number of iterations in case we never reach the break case.
t = ones(N,1); % Make a vector with N number of ones. The answers 
               % of every iteration of the newton method
               % will be stored in this vector.

x0 = 1; % Initial value assigned to x0. 
tol = 1e-5; % Same as 10^(-5). This is the tolerance level provided to us by the problem. 

% For loop that will go until N iterations or until the breaking case is reached.
for i = 1:N
   t(i) = x0 - f(x0)/fderi(x0); % Compute xn by using xn-1.
   if abs(t(i)-x0) <= tol % Check the stoping criteria/tolerance level.
      break; % Break if the tolerence level was met.
   end
   x0 = t(i); % Update the variable xn-1 to xn.
end

disp( 'Here are the 4 iterations for the Newton Method. NOTE: x0 is NOT shown since it is not part of the iteration process.' ) 
% Printing the list of numbers that got closer to the real root. 
t(1:4)
disp( '' ) % Print a newline. 

% The total number of iterations was 4. NOTE: x0 does not count as an
% iteration.
disp( 'How many iterations was run?' ) 
disp( '(Note: We do not count x0 as an iteration)' )
disp( 'The total number of iterations was the following:' ) 
i

% (c) Use the Matlab command r = fzero(fun, 1) to find the actual root, where fun is the
% function handle that defines 2^(-x) - x. Then print the first 5 (0, 1, 2, 3, 4) absolute errors
% for both methods. Use the long format to show 15 digits after the decimal point.

format long; % Make the format long to display 15 digits after the decimal point.
fun = @(x) 2^(-x) - x; % Define the new function fun.
r = fzero(fun, 1); % Find the actual root.

% Create two vectors that will store the 5 absolute errors of both methods.
pError = ones(5, 1);
tError = ones(5, 1);

% Since p0 = .5 and x0 = 1, manually compute the n=0 errors. This is
% required to be manual since x0 and p0 are not stored in the original
% iteration vectors.
pError(1) = abs( .50 - r );
tError(1) = abs( 1 - r );

% Use a for loop that will compute the relative errors of both methods from
% n = 1 to n n = 4.
for i = 2:5
   
    % Compute absolute errors. NOTE: We must use i-1 to access the
    % iteration roots since these vectors start at 1, not 2.
   pError(i) = abs( p(i -1) - r );
   tError(i) = abs( t(i -1) - r );
   
end

% Combine both absolute errors into a single matrix.
twoColumnErrors = [pError, tError];

% Print the errors side to side with a meaningful message.
disp( 'First five absolute errors. First column are the absolute errors of the Bisection Method.' )
disp( 'Second column are the absolute errors of Newtons method. Newtons Method produces better approximations faster.' ) 
twoColumnErrors


%% Ch2-2

% Newton's method fails for finding a root of y = x^5 - x - 1 with x0 = 0. Print the first 6
% approximations to illustrate this.

y = @(x) (x^5) - x - 1; % Define a new function for y.
yderi = @(x) 5*(x^4) - 1; % Declare the FIRST derivate of function f.

N = 6; % Set the max number of iterations to 6.
z = ones(N,1); % Make a vector with N number of ones. The answers 
               % of every iteration of the Newton method
               % will be stored in this vector.

x0 = 0; % Initial value assigned to x0. Note: This DOES NOT count as an iteration. 

% For loop that will run for 6 iterations. 
for i = 1:N
   z(i) = x0 - y(x0)/yderi(x0); % Compute xn by using xn-1.
   x0 = z(i); % Update the variable xn-1 to xn.
end

% Display the 6 iterations along with a meaningful message.
disp( 'Here are the first 6 iterations of Newtons Method for the given function.') 
disp( 'As the iterations demonstrate, Newtons Method fails since there is no convergence.')
disp( 'In fact, the iterations illustrate an oscillation rather than a convergence.')

z
