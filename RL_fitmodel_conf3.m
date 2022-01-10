function [loglik, pout]= RL_fitmodel_conf3(params, data, choice)
% FUNCTION [loglik pout]= RLtutorial_fitmodel(params, data)
%
% params = vector of parameter values for alpha and beta
% data.choice = vector of choices [1 2]
% data.outcome = vector of outcomes received, [0 1]
%
% The learning model used to fit choices is a simple Rescorla-Wagner
% (Rescorla & Wagner 1972) function, using a softmax
% choice/link/observation function:
%   learning model (RW):     vA <- vA + alpha*(r-vA)
%   choice function (softmax):    p(A) = exp(beta*vA)/(exp(beta*vA)_exp(beta*vB))
%
% ------------------------------------------------------------------------
% Adapted from a script written by Hanneke den Ouden in 2015 <h.denouden@gmail.com>
% Donders Center for Cognitive Neuroimaging
% Donders Center for Brain, Cognition and Behavior
% Radboud University Nijmegen
% Author: Elia Benhamou, University College London, UK - 2019
% ------------------------------------------------------------------------


% get behavioural data
outcome = data.outcome (data.prep.conf==3);

nt      = length(choice); %20 trials for each conf
ns      = 2; % number of choice options

% get params
alpha   = params(1);% learning rate
beta    = params(2); % sensitivity to reward and punishment
R       = params(3); %reward sensitivity 

% initialise variables
v0            = .5*ones(1,2); % initial value

r0            = 1;              %initial reward magnitude 
VVs          = nan(nt,ns); % matrix to store value over trials ns: for the 3 conf
PPs          = nan(size(VVs)); % matrix to store choice probability over trials

pe           = nan (nt, 1); 


% get values from all trials

v = v0;

for t = 1:nt
    c  = choice(t);
    r   = outcome(t);    
    % compute likelihood of the observed choices given the params and the
    % probability of 'go'

    ev  = exp(beta*v); % expected value
    sev = sum(ev);
    p   = ev/sev;
    
    % store value of V and choice probability

    VVs(t,:) = v; % store value of V
    PPs(t,:) = p;    % update value
    
    
    dv  = R*r-v(c+1); % compute prediction error
    v(c+1)   = v(c+1) + alpha*dv;    % update value
    
    %save prediction errors
    pe (t) = dv; 

    
end

loglik = sum(log(PPs(choice==0,1)))+sum(log(PPs(choice==1,2)));

% output parameters
pout = struct('PP3s',PPs, 'VV3s', VVs, 'loglik3',loglik, ...
     'pe3', pe,  'params3',params,'data',data);
end


