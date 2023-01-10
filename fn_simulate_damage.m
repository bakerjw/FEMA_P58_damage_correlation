function [results] = fn_simulate_damage(nSims, edp, component)
% Demo implementation of new component correlation algorithm
% created by Jack Baker, 4/20/2021
% updated 7/5/2022 to make paper-formatted composite figure
%
% INPUTS
%  nSims - how many Monte Carlo simulations are desired
%  edp - structure with fields describing the edp distribution
%  component - structure with characteristics and quantities of the components--a simple version of the real structure
%  makeFigs - plot figures if this = 1
%
% OUTPUTS
%  results - structure with some example outputs (to facilitate comparitive runs)
%

% make a/b/c labels for subplots
alphabet = ('a':'z').';
chars = num2cell(alphabet);
chars = chars.';
charlbl = strcat('(',chars,')'); % {'(a)','(b)','(c)','(d)'}
labelSize = 9;


%% weight cases to consider
w = [0   0  1   ; ...
     0.5 0  0.5 ; ...
     1   0  0  ];

nCases = size(w,1);


%% simulate EDPs (re-use these for all cases)
muEDP = log(edp.theta*ones(edp.nEDP,1)); % mean log EDPs
sigmaEDP = edp.beta*ones(edp.nEDP,1); % log standard deviations of EDPs

RHO = edp.rhoEDP*ones(edp.nEDP, edp.nEDP) + (1-edp.rhoEDP)*eye(edp.nEDP); % EDP correlation matrix
SIGMA = diag(sigmaEDP) * RHO * diag(sigmaEDP); % covariance matrix

edpSims = exp(mvnrnd(muEDP, SIGMA, nSims)); % realizations from target distribution (lognormal with muEDP and SIGMA)



%% simulate component variates
figure % start the figure, and add subplots inside the loop

for k = 1:nCases % for each weight case
    
    weights = w(k,:) / sum(w(k,:)); % pull the kth row of the weight matrix, and renormalize if weights didn't sum to 1

    variate_Total = normrnd(0,1,nSims,1); % one variate per sim that is used for all components
    variate_CompType = normrnd(0,1,nSims, length(component)); % one variate per component type

    nDamaged = zeros(nSims, length(component)); % initialize matrix
    for i = 1:length(component) % for each component type


        variate_i = normrnd(0, 1,nSims,edp.nEDP); % one component variate per EDP (assume all components with a given EDP are perfectly correlated, for now)

        % compute the component capacity
        compVariate = sqrt(weights(1))*repmat(variate_Total,1,edp.nEDP) + sqrt(weights(2))*repmat(variate_CompType(:,1),1,edp.nEDP) + sqrt(weights(3))*variate_i;
        compCapacities = exp(log(component(i).theta) +  component(i).beta * compVariate);

        % check whether each component is damaged, and count up numbers of damaged components
        isDamaged = edpSims > compCapacities;
        temp = isDamaged .* repmat(component(i).quantity, nSims, 1);
        nDamaged(:,i) = sum(isDamaged .* repmat(component(i).quantity, nSims, 1),2);
    end


    nDamagedTotal = sum(nDamaged, 2);

    % numerical outputs
    results.pDamage(k) = sum(nDamagedTotal>0)/nSims;
    results.p10(k) = sum(nDamagedTotal>9)/nSims;
    results.nDamaged(k) = mean(nDamagedTotal);


    % num damaged components histogram
    subplot(nCases, 2, k*2-1)
    histogram(nDamagedTotal);
    axis([-0.5 20.5 0 nSims])
    xlabel('Number of damaged components (N)', 'Interpreter','latex')
    ylabel('Number of simulations')
    %     title(['fractions = [' num2str(c(1),2) ',' num2str(c(2),2) ',' num2str(c(3),2) '], \mu_n=' num2str(mean(nDamagedTotal),2) ', \sigma_n=' num2str(std(nDamagedTotal),2)])
    axis square
    text(0.85,0.9,charlbl{k*2-1},'Units','normalized','FontSize',labelSize)
    text(0.15,0.8,['$P(N \geq 1) = ' num2str(results.pDamage(k),2) '$'],'Units','normalized','FontSize',labelSize, 'Interpreter','latex')
    text(0.15,0.7,['$P(N \geq 10) = ' num2str(results.p10(k),2) '$'],'Units','normalized','FontSize',labelSize, 'Interpreter','latex')
    text(0.15,0.6,['$E[N] = ' num2str(results.nDamaged(k),2) '$'],'Units','normalized','FontSize',labelSize, 'Interpreter','latex')


    % num damaged components scatter plot
    nP = 1000; % how many points to show

    jitters = (rand(nSims,2)-.5)*0.15; % random offsets for visualization

    subplot(nCases, 2, k*2)
    plot(nDamaged(1:nP,1)+jitters(1:nP,1), nDamaged(1:nP,2)+jitters(1:nP,2), '.')
    xticks(0:5)
    yticks(0:5)
    xlabel('Num. damaged components (type 1)')
    ylabel('Num. damaged components (type 2)')
    axis square
    axis([[0 length(component)+1 0 length(component)+1]-0.5])
    %     title(['\rho = ' num2str(results.damageRho,2)])
    text(0.85,0.9,charlbl{k*2},'Units','normalized','FontSize',labelSize)

end


% Format the figure sizes and fonts

% target font sizes
labelSize = 9;
axisSize = 9;

% use fixed figure size
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6 8]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 6 8]);

% loop over any subplots that exist within
panels = get(gcf,'Children');
for p = 1:length(panels)

    % check on type of child attribute
    if strcmp(get(panels(p),'type'), 'axes')

        % resize axis numbers
        set(panels(p), 'FontSize', axisSize);

        % resize axis labels
        axLabels = get(panels(p),{'XLabel', 'YLabel', 'ZLabel'});
        set([axLabels{:}], 'FontSize', labelSize);
    end
end
% FormatSubplotFigureBook
print(['simple_example.pdf'],'-dpdf')



end




