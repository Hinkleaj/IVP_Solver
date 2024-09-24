%==========================================================================
% PROJECT FOR NMoS: INVESTIGATION OF THE INFLUENCE OF AN ADJUSTABLE 
%                   STEP-SIZE IN THE BACKWARD EULER METHOD
%--------------------------------------------------------------------------
%Three methods of adjusting h are compared in their computational time
%accuracy and stability:
%   1. constant h 
%   2. using the second derivative
%   3. adjusting the h in a way to minimize the MSE by using gradient
%       descent
%
%To test these methods functions with a known solution are tested.
%--------------------------------------------------------------------------
%In this file the Implicit Euler MEthod (IEM) is done with adjusting the h 
% in a way to minimize the MSE by using gradient descent
%==========================================================================


%--------------------------------------------------------------------------
% PUT IN INITIAL VALUES, FUNCTION AND BOUNDARIES
%--------------------------------------------------------------------------
clear all;
x_start = 5;
t_start = 0;
t_end = 5;
h_start = 0.1;
step_newt = 10000;
error = 0.001;

%--------------------------------------------------------------------------
% START SOLVING WITH ADJUSTMENT BY MINIMIZING MSE
%--------------------------------------------------------------------------
t_array(1) = t_start;
i = 1;
h(1) = h_start;
while t_array(i)<=t_end %LOOP FOR EULER STEPS
   x_sol(1) = x_start;
   
   h_new = gradient_descent(@MSE,h_start,x_sol(i),t_array(i),@dx_dt_sol,x_start,step_newt);
   h(i+1)=h_new;
   x_sol(i+1)=newton(@dx_dt,x_sol(i),x_sol(i)+0.001,t_array(i)+h_new,h_new,step_newt);
   t_array(i+1)=t_array(i)+h_new;
   i = i+1;

end %LOOP FOR EULER STEPS

%--------------------------------------------------------------------------
% REAL VALUES OF X AND MSE
%--------------------------------------------------------------------------

for j=1:i-1 %LOOP FOR REAL VALUES
   x_true(1) = x_start; %solution of x with the newton method 
   mse_value(1) = 0;
   x_true(j+1) = dx_dt_sol(t_array(j+1),x_start);
   mse_value(j+1) = MSE(x_sol(j+1),x_true(j+1));
end %LOOP FOR REAL VALUE
mse_sum = sum(mse_value);
%--------------------------------------------------------------------------
% Visualisation of the Solution
%--------------------------------------------------------------------------
tiledlayout(1,2)
nexttile
plot(t_array,x_sol,"--",'Color',[0, 0.5, 0.3],'LineWidth',1);
hold on
plot(t_array,x_true,"--",'Color',[0.5, 0, 0.3],'LineWidth',1);
xlabel("t");
ylabel("x");
nexttile
plot(t_array,mse_value,"--",'Color',[0.3, 0, 0.5],'LineWidth',1);





%==========================================================================
% FUNCTIONS
%==========================================================================

%--------------------------------------------------------------------------
% FUNCTION TO SOLVE AND FUNCTION FOR REAL SOLUTION (CHANGE THESE FOR YOUR
% PuRPOSE)
%--------------------------------------------------------------------------
function func = dx_dt(t,x)
  func = 2*x;
end

function f_solution = dx_dt_sol(t,x_0)
    f_solution = x_0*exp(2*t);
end

%--------------------------------------------------------------------------
% FUNCTION FOR GRADIENT DESCENT
%------------------------------------------------------------------------

function h_right = gradient_descent(funct,h_start,x_sol,t,x_true,x_start,step_newt)
    min(1)= h_start;
    k = 2;
    alpha = 0.01;
    f_prime = (funct(newton(@dx_dt,x_sol,x_sol+0.001,t+min(1)+0.00001,min(1)+0.00001,step_newt),x_true(t+min(1),x_start))-funct(newton(@dx_dt,x_sol,x_sol+0.001,t+min(1),min(1),step_newt),x_true(t+min(1),x_start)))/0.00001;
    min(k) = min(1) - alpha * f_prime;
    while abs(min(k-1)-min(k))>0.1 & k < 10000
        f_prime = (funct(newton(@dx_dt,x_sol,x_sol+0.001,t+min(k)+0.00001,min(k)+0.00001,step_newt),x_true(t+min(k),x_start))-funct(newton(@dx_dt,x_sol,x_sol+0.001,t+min(k),min(k),step_newt),x_true(t+min(k),x_start)))/0.00001;
        min(k+1) = min(k) - alpha * f_prime;
        k=k+1;
    end
    if min(k)>0 & abs(min(k-1)-min(k))<0.1
        h_right=min(k);
    else
        h_right=0.001;
    end
    
end

%--------------------------------------------------------------------------
% FUNCTION FOR NEWTON METHOD
%--------------------------------------------------------------------------
function x_newt = newton(dx_dt,x_i, x_guess, t, h,step_newt)
  prime_factor = 0.00000001; %For getting the drivative
  error = 0.0001; %Tolerance for finding the root
  %First newton step
  x_search(1) = x_guess;
  
  funct = dx_dt(t,x_guess);
  newton =  x_i + h * funct - x_guess;
  primestep = (x_i + h * dx_dt(t,x_guess+prime_factor) - (x_guess+prime_factor));
  prime = (primestep-(newton))/prime_factor;
  j = 2;
  x_search(j) = x_search(1)-newton/prime;
  while (abs(x_search(j)-x_search(j-1))>error) & j < step_newt %LOOP FOR NEWTON STEPS
    funct = dx_dt(t,x_search(j));
    newton =  x_search(j) + h * funct - x_search(j);
    prime = ((x_i + h * dx_dt(t,x_search(j)+prime_factor) - (x_search(j)+prime_factor))-(newton))/prime_factor;
    x_search(j+1)=x_search(j)-newton/prime;
    j=j+1;
  end%LOOP FOR NEWTON STEPS
  if abs(x_search(j)-x_search(j-1))>error%IF STATEMENT FOR NEWTONSTEPS (if not successful)
      x_newt  = 0;
  end%IF STATEMENT FOR NEWTONSTEPS
  x_newt  = x_search(j); 
end

%--------------------------------------------------------------------------
% FUNCTION FOR MEAN SQUARE ERROR
%--------------------------------------------------------------------------

function MSE = MSE(x_estimate,x_true)
    MSE = (x_estimate-x_true)^2;
end