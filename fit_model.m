% function fit_model(expid, subjidx)
%
% This function fits the resource-rational model as described in the paper:
%
%   *********************************************************************
%   * Van den Berg, R. & Ma, W.J. (2018). A resource-rational theory of *
%   *   set size effects in human visual working memory. Elife.         *
%   *********************************************************************
%
% INPUT
%  expid   : experiment ID (integer in range 1 to 9; see Table 1 in paper)
%  subjidx : subject index (integer in range 1 to #subjects)
%
% Written by Ronald van den Berg, 2018

function fit_model(expid, subjidx)

if expid==8 || expid==9
    fprintf('\n---------------------------------------------------------------\n');
    fprintf('Sorry, due to legal restrictions we were not allowed to publish\n')
    fprintf('the datasets from the study by Emrich et al. here. The data can\n');
    fprintf('be requested by contacting Stephen Emrich at semrich@brocku.ca.\n');
    fprintf('---------------------------------------------------------------\n');
    return
end

if ~exist('bads','file')
    fprintf('\n');
    fprintf('----------------------------------------------------------------------------------------\n');
    fprintf('This function requires the Bayesian Adaptive Direct Search (BADS) optimization algorithm\n');
    fprintf('for model fitting in MATLAB. Please visit https://github.com/lacerbi/bads and press the \n');
    fprintf('green button that says "Clone or Download" to get a copy.                               \n');
    fprintf('                                                                                        \n');
    fprintf('Alternatively, you could use fminsearch, by replacing the BADS call (line 66) with:     \n');
    fprintf('   [fitpars, mLLH] = fminsearch(@(pars) -LLH_function(pars, data, gvar), fp_init)       \n');
    fprintf('                                                                                        \n');
    fprintf('However, BADS is recommended as it is faster and less likely to end in a local minimum  \n');
    fprintf('----------------------------------------------------------------------------------------\n');
    return
end

% Load data
load(sprintf('data/data_E%d.mat',expid),'data');
data = data{subjidx};

% Set up "global" variable structure for easy passing to functions
gvar.kmap = [linspace(0,10,250) linspace(10.001,20000,500)];
gvar.Jmap = gvar.kmap.*besseli(1,gvar.kmap,1)./besseli(0,gvar.kmap,1); % mapping from kappa to J (see Appendix 1)
gvar.n_gamma_bins = 50;        % number of bins to use to discretize Gamma distribution over J when computing model predictions
gvar.n_VM_bins = 360;          % number of bins to use to discretize Von Mises distribution over estimation error when computing model predictions
gvar.n_initpars_try = 25;      % number of random parameter vectors to try for initialization of optimizer
gvar.uPi = unique(data.p_i)';  % unique values of p_i in this dataset

fprintf('\nNOTE: Increase n_gamma_bins, n_VM_bins, and n_initpars_try to get more precise results\n\n');

% Find initial parameter estimate for optimizer
fprintf('Searching for initial parameter setting...\n');
[fp_init, lb, ub, plb, pub, maxLLH_init] = get_initpars(data,gvar);
fprintf(' Maximum log likelihood (init) = %2.1f\n',maxLLH_init);
fprintf(' Parameter estimates    (init) = %2.2f (tau), %2.5f (lambda), %2.3f (beta)\n\n',fp_init(1),fp_init(2),fp_init(3));

% Run optimizer to find maximum-likelihood parameter estimates
fprintf('Running optimizer to find maximum-likelihood parameter estimates...\n');
bads_options = bads('defaults');
bads_options.Display = 'off';
[fitpars, mLLH] = bads(@(pars) -LLH_function(pars, data, gvar), fp_init, lb, ub, plb, pub, bads_options);
fprintf(' Maximum log likelihood (final) = %2.1f\n',-mLLH);
fprintf(' AIC                    (final) = %2.1f\n',2*mLLH + 2*numel(fitpars));
fprintf(' Parameter estimates    (final) = %2.2f (tau), %2.5f (lambda), %2.3f (beta)\n',fitpars(1),fitpars(2),fitpars(3));

% plot result
plot_data_and_fit(data,fitpars,gvar);

%------------------- HELPER FUNCTIONS -------------------%

% This function returns the log likelihood for a given set of parameters and dataset
function LLH = LLH_function(pars,data,gvar)
tau=pars(1);
lambda=pars(2);
beta=pars(3);
uPi = unique(data.p_i);
p_resp = zeros(1,numel(data.p_i));
for ii=1:numel(uPi)
    % Find indices of all trials with this value of p_i
    idx = find(data.p_i==uPi(ii));            
    % Compute optimal Jbar for this value of p_i
    Jbar_optimal = fminsearch(@(pars) cost_function(pars,tau,lambda,beta,uPi(ii),gvar), 1);
    % Compute probability of the subject's estimation errors under this value of Jbar_optimal
    J = discretize_gamma(Jbar_optimal,tau,gvar.n_gamma_bins);
    kappa = interp1(gvar.Jmap,gvar.kmap,min(J,max(gvar.Jmap)));
    p_resp(idx) = mean(bsxfun(@times,1./(2*pi*besseli(0,kappa)),exp(bsxfun(@times,kappa,cos(data.error(idx))))),2);
end
LLH = sum(log(max(p_resp,1e-3)));

% This function returns expected total cost for a given set of parameters and p_i value
function E_C_total=cost_function(Jbar,tau,lambda,beta,p_i,gvar)
J = discretize_gamma(Jbar,tau,gvar.n_gamma_bins);
kappa = interp1(gvar.Jmap,gvar.kmap,min(J,max(gvar.Jmap)));
VM_x = linspace(-pi,pi,gvar.n_VM_bins);
VM_x = VM_x(2:end)-diff(VM_x(1:2))/2;
VM_y = exp(bsxfun(@times,kappa',cos(VM_x)));
VM_y = bsxfun(@rdivide,VM_y,sum(VM_y,2));
C_behavioral = abs(VM_x).^beta; 
E_C_behavioral = mean(sum(bsxfun(@times,VM_y,C_behavioral),2)); 
E_C_total = p_i*E_C_behavioral + lambda*Jbar; 

function [initpars, lb, ub, plb, pub, maxLLH] = get_initpars(data,gvar)
% parameters: [tau, lambda, beta]
lb  = [1e-4  1e-8    0];   % lower bounds on parameters
ub  = [ 200     1   10];   % upper bounds on parameters 
plb = [  .1  1e-6    0];   % plausible lower bounds
pub = [  10   .01    3];   % plausible upper bounds
% try a bunch of randomly picked parameter values and return the best one
maxLLH = -Inf;
for ii=1:gvar.n_initpars_try 
    parvec = plb+rand*(pub-plb);
    LLH = LLH_function(parvec,data,gvar);
    if LLH>maxLLH
        maxLLH = LLH;
        initpars = parvec;
    end
end

% This function creates a plot of the data together with model fit
function plot_data_and_fit(data,fitpars,gvar)
figure
set(gcf,'Position',get(gcf,'Position').*[0.1 0.1 1.5 1]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[0.1 0.1 1.5 1])
tau=fitpars(1);
lambda=fitpars(2);
beta=fitpars(3);
uPi = unique(data.p_i);
nbins_emp = 25;
nbins_fit = 199;
X_emp = linspace(-pi,pi,nbins_emp+1);
X_emp = X_emp(2:end)-diff(X_emp(1:2))/2;
X_fit = linspace(-pi,pi,nbins_fit+1);
X_fit = X_fit(2:end)-diff(X_fit(1:2))/2;
plotids = zeros(1,8);
for ii=1:numel(uPi)
    % compute empirical error distribution for this p_i value
    Y_emp = hist(data.error(data.p_i==uPi(ii)),X_emp);
    Y_emp = Y_emp/sum(Y_emp)/diff(X_emp(1:2));    
    % compute error distribution corresponding to model fit
    Jbar = fminsearch(@(pars) cost_function(pars,tau,lambda,beta,uPi(ii),gvar), 1);
    J = discretize_gamma(Jbar,tau,gvar.n_gamma_bins);
    kappa = interp1(gvar.Jmap,gvar.kmap,min(J,max(gvar.Jmap)));
    Y_fit = mean(bsxfun(@times,1./(2*pi*besseli(0,kappa)),exp(bsxfun(@times,kappa,cos(X_fit')))),2);        
    % plot
    plotids(ii) = subplot(2,4,ii);
    hold on
    plot(X_emp,Y_emp,'ko','markerfacecolor','k','markersize',4);
    plot(X_fit,Y_fit,'k-');
    title(sprintf('p_i=%2.2f',uPi(ii)));
end
warning off; linkaxes(plotids,'y'); warning on;

% This function discretizes a gamma distribution and returns bin centers
function bins = discretize_gamma(Jbar,tau,nbins)
X = linspace(0,1,nbins+1);
X = X(2:end)-diff(X(1:2))/2;
warning off
bins = gaminv(X,Jbar/tau,tau);
warning on

