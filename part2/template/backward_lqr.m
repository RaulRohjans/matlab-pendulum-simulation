function [d,K] = backward_lqr(X,U,A,B,Q,R,Qf)
% A is of size (4,4,N), A matrices corresponding to linearized trajectory
% B is of size (4,1,N), B matrices correspondng to linearized trajectory


N = numel(U);  %number of time steps

%%% Initialization
d = zeros(1,N); % feedforward terms 
% Note: d corresponds to the feedforward k term in the class notes

K = zeros(1,4,N); %state feedback gains
 
%%% Implement backward pass here 

end