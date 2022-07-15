function [c, ceq, DC, DCeq] =  fuel_beta_constraint(variables)

% This is + keff because the python code outputs it as a negative value
c = sum(1 - variables) - 57;
ceq = 0;
DC = ones(121,1) * -1;
DCeq = zeros(121,1);

gamma = ones(121, 1) * 1;

penalty = gamma*sqrt(sum(variables.*(1-variables)))
penalty_grad = gamma*(-2*variables+ones(size(variables)))/(2*sqrt(sum(variables.*(1-variables))))
variables
penalty(1)

c1 = sum(1 - variables) - 57;
c2 = penalty(1)
c = [c1; c2]
% The penalty deriv is added to the original deriv, -1


DC1 = ones(121,1) * -1
DC2  = penalty_grad(1,:);
DC2  = reshape(DC2, 121, 1)
DC = [DC1; DC2]
DC = reshape(DC, 121, 2)