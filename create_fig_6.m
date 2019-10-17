% function create_fig_3A_bottom
%
% This function creates Figure 6 of the paper:
%
%   *********************************************************************
%   * Van den Berg, R. & Ma, W.J. (2018). A resource-rational theory of *
%   *   set size effects in human visual working memory. Elife.         *
%   *********************************************************************
%
% Written by Ronald van den Berg, 2018

function create_fig_6

% Simulation settings
lambda = .0082;  % value of model parameter lambda
tau    = 30;     % value of model parameter tau
N      = 5;      % set size at which the experiment is performed

% Store "global" parameters in a structure variable for easy passing to functions
gvar.kmap = [linspace(0,10,250) linspace(10.001,20000,500)];
gvar.Jmap = gvar.kmap.*besseli(1,gvar.kmap,1)./besseli(0,gvar.kmap,1); % mapping between J and kappa (see Appendix 1)
gvar.n_gamma_bins = 50;  % number of bins to use for discretization of gamma distribution over J when computing model predictions
gvar.n_VM_bins = 720;    % number of bins to use for discretization of Von Mises distribution over estimation error when computiong model predictions
gvar.nsamples = 1e5;     % number of MC samples to use for model predictions (INCREASE THIS TO GET SMOOTHER PLOT)
gvar.nthresholds = 50;   % number of threshold values for which to comptue predictions (INCREASE THIS TO GET SMOOTHER PLOT)

fprintf('\nTo get a smoother plot, increase variables gvar.nsamples and gvar.nthresholds\n\n');

% Compute model predictions
p_i = 1/N;
th_vec = linspace(0,pi,gvar.nthresholds); % vector with "reward thresholds" for which to compute model predictions
Y_mean = zeros(1,numel(th_vec));
Y_std = zeros(1,numel(th_vec));
fprintf('0%% %s 100%%\n',blanks(gvar.nthresholds-8));
for ii=1:numel(th_vec)
    fprintf('.');
    Jbar_opt = fminsearch(@(Jbar) cost_function(Jbar,tau,lambda,th_vec(ii),p_i,gvar), 2);
    % computed mean(abs(error)) under this Jbar
    J = gamrnd(Jbar_opt/tau,tau,1,gvar.nsamples);
    kappa = interp1(gvar.Jmap,gvar.kmap,J,'linear','extrap');
    x = circ_vmrnd(0,kappa);
    % compute mean absolute error for plotting 
    Y_mean(ii) = mean(abs(x));
    Y_std(ii) = std(abs(x));
end
fprintf('\n')

% Rest is plotting
figure
hold on
plot([pi/8 pi/3],[1.025 .96],'ko','markerfacecolor','k','markersize',8);
plot(th_vec,Y_mean,'r-');
set(gca,'Xtick',[0 pi/4 pi/2 3/4*pi pi],'Xticklabel',{'0','\pi/8','\pi/2','3\pi/4','pi'},'Ytick',[0 pi/6 pi/3 pi/2],'YtickLabel',{'0','pi/6','pi/3','pi/2'});
plot([pi/3 pi/3],[0 1.6],'k--');
plot([pi/8 pi/8],[0 1.6],'k--');
xlim([0 pi]);
ylim([0 pi/1.9])
ylabel('Average absolute estimation error');
xlabel('Positive feedback threshold');
grid on


% The below function returns expected total cost 
function E_C_total=cost_function(Jbar,tau,lambda,th,p_i,gvar)
J = discretize_gamma(Jbar,tau,gvar.n_gamma_bins);  % discretize 
kappa = interp1(gvar.Jmap,gvar.kmap,min(J,max(gvar.Jmap)));
VM_x = linspace(-pi,pi,gvar.n_VM_bins);
VM_x = VM_x(2:end)-diff(VM_x(1:2))/2;
VM_y = exp(bsxfun(@times,kappa',cos(VM_x)));
VM_y = bsxfun(@rdivide,VM_y,sum(VM_y,2));
C_behavioral(abs(VM_x)<th)=0;
C_behavioral(abs(VM_x)>th)=1;
E_C_behavioral = mean(sum(bsxfun(@times,VM_y,C_behavioral),2)); % expected behavioral cost
E_C_total = p_i*E_C_behavioral + lambda*Jbar; % local expected total cost

function bins = discretize_gamma(Jbar,tau,nbins)
X = linspace(0,1,nbins+1);
X = X(2:end)-diff(X(1:2))/2;
warning off
bins = gaminv(X,Jbar/tau,tau);
warning on

% This is a modified copy of the circ_vmrnd function from the Circular Statistics
% Toolbox for Matlab (by Philipp Berens and Marc J. Velasco, 2009)
function Y = circ_vmrnd(theta, kappa)
n=1;
input_dims = size(kappa);
kappa = kappa(:)';
kappa = repmat(kappa,1,n);
theta = ones(size(kappa))*theta;
a = 1 + sqrt((1+4*kappa.^2));
b = (a - sqrt(2*a))./(2*kappa);
r = (1 + b.^2)./(2*b);
valid = zeros(1,length(kappa));
z = zeros(size(kappa));
f = zeros(size(kappa));
c = zeros(size(kappa));
while ~all(valid)
    u(:,~valid) = rand(3,sum(~valid));    
    z(~valid) = cos(pi*u(1,~valid));
    f(~valid) = (1+r(~valid).*z(~valid))./(r(~valid)+z(~valid));
    c(~valid) = kappa(~valid).*(r(~valid)-f(~valid));       
    valid = u(2,:) < c .* (2-c) | ~(log(c)-log(u(2,:)) + 1 -c < 0);               
end
Y = theta + sign(u(3,:) - 0.5) .* acos(f);
Y = angle(exp(1i*Y));
Y(kappa<1e-6) = 2*pi*rand(sum(kappa<1e-6),1)-pi;
Y = reshape(Y,input_dims);

