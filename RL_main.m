function [score_levels, RT_levels] = RLtutorial_main_diffconf
% Adapted from a script written by Hanneke den Ouden in 2015 <h.denouden@gmail.com>
% Donders Center for Cognitive Neuroimaging
% Donders Center for Brain, Cognition and Behavior
% Radboud University Nijmegen
% Author: Elia Benhamou, University College London, UK - 2019

%% Section 1: Preparation
%==========================================================================
close all;

%%%%%%%%%%%%    MODIFY      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use vector of numbers if simulating or using a subset of the data; use string 'all' if using all available real datasets
subjects        = [11 21 31 41 51 61 71 81 91 111 121 131 141 151 161 171 191 201 211 221 231];


% CONS controls: [11 21 31 41 51 71 81 91 101 111 121 131 141 161 171 181  201 211 221 231 241 251 261] 
% SQ controls: [13 23 33 43 53 63 83 93 103 113 123 133 143 153 163 173 183 193];
% CONS bv: [11 21 31 41 51 61 71 81 91 111 121 131 141 151 161 171 191 201 211 221 231]; pikge: 61
% SQ bv: [13 23 33 43 53 63 83 93 103 113 123 133 143 163 173 183];
% CONS ad: [11 31 41 51 61 71 81 91 101 111];
% SQ ad: [23 33 43 53 63 73 83 93];

simulate        = false;    % true if simulate, false if using real data
fitData         = false; 	% if false, we only simulate / plot data
plotIndividual  = false;  % only set to to true if running a small number of individuals!

% If simulating data, specify the parameters you want to use for fitting
if simulate
    simpars.alpha = 0.25;
    simpars.beta  = 4;
end

% If fitting data, specify the bounds of the grid, and the number of bins
% that we want to compute.
if fitData
    bounds1      = [0 1; % alpha
        0 25;  % beta
        0 1]; % R
    bounds2      = [0 1; % alpha
        0 25;  % beta
        0 1]; % R
    bounds3      = [0 1; % alpha
        0 25;  % beta
        0 2]; % R
    nBin        = [40 50 50] ;
end

%%%%%%%%%%%%    MODIFY      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set paths for code we need
rootdir = '/Users/Eliacello/Documents/PhD_DRC/RLmodel/RLtest1/';
addpath(fullfile(rootdir,'util')); % some helper functions
addpath(fullfile(rootdir,'data'));% to get the feedback

% code-rewrite
dataPath = fullfile(rootdir,'data'); % where our data comes from
dataPath_bv = fullfile(rootdir,'data_bv');
dataPath_ad = fullfile(rootdir,'data_ad');

% set labels that we'll need later on.
level  = {'90% vs 10%','90% vs 50%','50% vs 10%'};
RT = {'90%','50%','10%'};
paramLabel  = {'alpha','beta', 'R'};
nParam      = length(paramLabel);


%==================================================================
%% Section 2: Simulate data or plot real data
%==================================================================

% initialise variables
score_all = nan (1, length(subjects));
score_easy = nan (1, length(subjects));
score_hard = nan (1, length(subjects));
RT_rew = nan (1, length(subjects));
RT_amb = nan (1, length(subjects));
RT_pun = nan (1, length(subjects));

% score           = nan(2,nsubTask);
if ~plotIndividual; fhData = figure('Name','Data');fhData1 = figure('Name','Data'); fhData2 = figure('Name','Data');
    set(fhData,'position',[10 60 900 650],'paperunits','centimeters','Color','w');
    set(fhData1,'position',[10 60 900 650],'paperunits','centimeters','Color','w');
    set(fhData2,'position',[10 60 900 650],'paperunits','centimeters','Color','w');
    
else
    trialfh = nan(1,length(subjects)); % figure handle for individual trialwise plots
    trialfh1 = nan(1,length(subjects));
    trialfh2 = nan(1,length(subjects));
end

ct      = 0; % reset subject counter

for sID = subjects'
    for ct = 1:length(subjects)

            subTag      = sprintf('bv_S_%d_data.mat',sID(ct));% controls_S_%d_data.mat or bv_S_
            dataFile = fullfile(dataPath_bv, subTag); %dataPath or dataPath_bv 
            load(dataFile, 'data');
            alldata{ct} = data;
            
            if plotIndividual
                [trialfh(ct), trialfh1(ct), trialfh2(ct)] = revlPlotData(sID(ct),dataPath);
            end
            choice(ct, 1:data.prep.nt) =( 3-data.choice)/2;
            feedbackprob(1:data.prep.nt, ct) = data.prep.feedbackprob;
            data.prep.conf(1:data.prep.nt, ct) = data.prep.conf;
            
            choice1_all(ct, 1:20) = (3-data.choice(data.prep.conf(:,ct)==1))/2;
            choice1_all(choice1_all==0.5)=0;
            choice2_all(ct, 1:20) = (3-data.choice(data.prep.conf(:,ct)==2))/2;
            choice3_all(ct, 1:20) = (3-data.choice(data.prep.conf(:,ct)==3))/2;
            choice3_all(choice3_all==0.5)=1;
            
            score_all(ct) = data.totalAcc/data.prep.nt;
            score_easy(ct) = data.AccEasy/40;
            score_hard(ct) = data.AccHard/20;
            
            score_pair1(ct)=length(find(data.choice(data.prep.conf(:,ct)==1)==1))/20;
            score_pair2(ct)=length(find(data.choice(data.prep.conf(:,ct)==2)==1))/20;
            score_pair3(ct)=length(find(data.choice(data.prep.conf(:,ct)==3)==2))/20;
            
            RT2=data.RT;
            RT2(data.RT>120000)=NaN;
            RT_rew(ct) = nanmedian(RT2(data.choice==1));
            RT_amb(ct) = nanmedian(RT2(data.choice==2));
            RT_pun(ct) = nanmedian(RT2(data.choice==3));
            
            RT_rewN(ct) = length(RT2(data.choice==1));
            RT_ambN(ct) = length(RT2(data.choice==2));
            RT_punN(ct) = length(RT2(data.choice==3));
    end
end



if ~plotIndividual
    figure(fhData)
    %     plot(data.prep.feedbackprob,'k:','linewidth',2)
    plot(nanmean(choice1_all,1),'k-','linewidth',1)
    ylabel('probability');
    ylim([0 1]);
    xlabel('trial');
    ylabel('probability')
    
    figure(fhData1)
    %     plot(data.prep.feedbackprob,'k:','linewidth',2)
    plot(nanmean(choice2_all,1),'k-','linewidth',1)
    ylabel('probability');
    ylim([0 1]);
    xlabel('trial');
    ylabel('probability')
    
    figure(fhData2)
    plot(nanmean(choice3_all,1),'k-','linewidth',1)
    ylabel('probability');
    ylim([0 1]);
    xlabel('trial');
    ylabel('probability')
    
 
        legend({'p(reward|all)','mean choice (rew)'},'location','northeastoutside');
        title(sprintf('MEAN data'));

    legend boxoff
    ylabel('probability');
    ylim([0 1]);
    xlabel('trial');
    ylabel('probability');
end % if plotting group

% plot the scores
if ~plotIndividual
    fhScore = figure('Name','score'); set(fhScore,'position',[10 60 250 400],'paperunits','centimeters',...
        'paperposition',[0 0 6 6],'Color','w');
    score_levels=[score_all; score_easy; score_hard];
    score_pairs=[score_pair2; score_pair1; score_pair3]; 
    barScatter(score_pairs,[],level,true, length(subjects));
    set(gca,'xtick',1:3,'xticklabel',level);
    xlabel('');
    ylabel('p(correct)');
    title('')
    hline(.5,'w:');
    ylim([0 1])
    box off
    saveas(gcf, 'scores_bv.pdf');
    
    fhRT = figure('Name','RT'); set(fhRT,'position',[10 60 250 400],'paperunits','centimeters',...
        'paperposition',[0 0 6 6],'Color','w');
    RT_levels=[RT_rew; RT_amb; RT_pun];
    barScatter(RT_levels,[],RT,true, length(subjects));
    set(gca,'xtick',1:3,'xticklabel',RT);
%     xlabel('arm');
    ylabel('Reaction time (ms)');
    %     ylim ([0 15000]);
%     title('reaction time')
    hline(.5,'w:');
    
    box off
    saveas(gcf, 'RT_bv.pdf');
    
    
end

%==========================================================================
%% Section 3.
% Do a grid search to compute the likelihood functions of the
% parameters, within prespecified bounds
%==========================================================================
if fitData
    if size(bounds1,1)~=nParam
        error('number of bounds and parameters don''t match')
    end
    
    % Initialise the output matrices for the fitted parameters.
    fitted.alpha.pdf    = nan(1,length(subjects),nBin(1));
    fitted.beta.pdf     = nan(1,length(subjects),nBin(2));
    fitted.R.pdf        = nan(1,length(subjects),nBin(3));
    fitted.alpha1.pdf    = nan(1,length(subjects),nBin(1));
    fitted.beta1.pdf     = nan(1,length(subjects),nBin(2));
    fitted.R1.pdf       = nan(1,length(subjects),nBin(3));
    fitted.alpha2.pdf    = nan(1,length(subjects),nBin(1));
    fitted.beta2.pdf     = nan(1,length(subjects),nBin(2));
    fitted.R2.pdf       = nan(1,length(subjects),nBin(3));
    
    fitted.alpha.ml     = nan(1,length(subjects));
    fitted.beta.ml      = nan(1,length(subjects));
    fitted.R.ml         = nan(1,length(subjects));
    fitted.alpha1.ml     = nan(1,length(subjects));
    fitted.beta1.ml      = nan(1,length(subjects));
    fitted.R1.ml        = nan(1,length(subjects));
    fitted.alpha2.ml     = nan(1,length(subjects));
    fitted.beta2.ml      = nan(1,length(subjects));
    fitted.R2.ml        = nan(1,length(subjects));
    
    fitted.alpha.ev     = nan(1,length(subjects));
    fitted.beta.ev      = nan(1,length(subjects));
    fitted.R.ev         = nan(1,length(subjects));
    fitted.alpha1.ev     = nan(1,length(subjects));
    fitted.beta1.ev      = nan(1,length(subjects));
    fitted.R1.ev        = nan(1,length(subjects));
    fitted.alpha2.ev     = nan(1,length(subjects));
    fitted.beta2.ev      = nan(1,length(subjects));
    fitted.R2.ev        = nan(1,length(subjects));
    
    % open some figures
    fhParamPDF = figure('Name','parameter PDF all');
    set(fhParamPDF,'position',[10 60 650 650],'paperunits','centimeters','Color','w');
    fhParamPDF1 = figure('Name','parameter PDF easy');
    set(fhParamPDF1,'position',[10 60 650 650],'paperunits','centimeters','Color','w');
    fhParamPDF2 = figure('Name','parameter PDF hard');
    set(fhParamPDF2,'position',[10 60 650 650],'paperunits','centimeters','Color','w');
    if plotIndividual
        fhGrid = figure('Name','Likelihoods  pair 1');
        set(fhGrid,'position',[10 60 650 650],'paperunits','centimeters','Color','w');
        fhGrid1 = figure('Name','Likelihoods pair 2');
        set(fhGrid1,'position',[10 60 650 650],'paperunits','centimeters','Color','w');
        fhGrid2 = figure('Name','Likelihoods pair 3');
        set(fhGrid2,'position',[10 60 650 650],'paperunits','centimeters','Color','w');
        
    else
        fhParBar = figure('Name','Parameter estimates');
        set(fhParBar,'position',[10 60 650 650],'paperunits','centimeters','Color','w');
        fhParBar1 = figure('Name','Parameter estimates');
        set(fhParBar1,'position',[10 60 650 650],'paperunits','centimeters','Color','w');
        fhParBar2 = figure('Name','Parameter estimates');
        set(fhParBar2,'position',[10 60 650 650],'paperunits','centimeters','Color','w');
    end
    
    %% Do the actual grid search
    % =========================================================================

    ct = 0; % subject counter
    
    for sID = subjects'
        for ct = 1:length(subjects)
            %                 ct = ct+1;
            %                 ctall = ctall+1;
            sID_string{ct} = sprintf('sub %03.0f',sID(ct));
            data = alldata{ct};
            choice1 = choice1_all(ct, :);
            choice2 = choice2_all(ct, :);
            choice3 = choice3_all(ct, :);
            
            for iParam = 1:nParam
                range1 = linspace(bounds1(iParam,1),bounds1(iParam,2),nBin(iParam)+1);
                p1{iParam} = range1(2:end); % stay just off the zero bounds

                range2 = linspace(bounds2(iParam,1),bounds2(iParam,2),nBin(iParam)+1);
                p2{iParam} = range2(2:end); % stay just off the zero bounds

                range3 = linspace(bounds3(iParam,1),bounds3(iParam,2),nBin(iParam)+1);
                p3{iParam} = range3(2:end); % stay just off the zero bounds

            end
            
            params1 = nan(1,3);
            params2 = nan(1,3);
            params3 = nan(1,3);
            for t = 1:nBin(1)
                params1(1) = p1{1}(t);
                params2(1) = p2{1}(t);
                params3(1) = p3{1}(t);
                for tt = 1:nBin(2)
                    params1(2) = p1{2}(tt);
                    params2(2) = p2{2}(tt);
                    params3(2) = p3{2}(tt);
                    for ttt = 1:nBin(3)
                        params1(3) = p1{3}(ttt);
                        params2(3) = p2{3}(ttt);
                        params3(3) = p3{3}(ttt);
                        
                    [loglik1(t,tt,ttt), foo]= RL_fitmodel_conf1(params1, data, choice1); %RL_fitmodel_3conf or RL_fitmodel_3conf
                    [loglik2(t,tt,ttt), foo]= RL_fitmodel_conf2(params2, data,  choice2);
                    [loglik3(t,tt,ttt), foo]= RL_fitmodel_conf3(params3, data, choice3);
                    end
                end
            end
            %             loglik1(loglik1==-Inf)=NaN;
            loglik1 = loglik1-min(loglik1(:)); % remove the minimum;
            %             loglik2(loglik2==-Inf)=NaN;
            loglik2 = loglik2-min(loglik2(:));
            %             loglik3(loglik3==-Inf)=NaN;
            loglik3 = loglik3-min(loglik3(:));
            lik1 = exp(loglik1); % compute the likelihood, (rather than the log)
            lik2 = exp(loglik2);
            lik3 = exp(loglik3);
            
            % compute the marginal likelihoods for each parameter: for
            % alpha, sum the likelihoods across all bins of beta and R,  for
            % beta, sum the likelihoods across all bins of alpha and R and for R across 
            % all bins of alpha and beta. Then to
            % compute the probability density, normalise by dividing by the sum
            % of the likelihoods for each parameter.
             x = 1; 
                tmp1 = sum(lik1,4-x); tmp1=sum(squeeze(tmp1),3-x);
                tmp2=sum(lik2,4-x); tmp2=sum(squeeze(tmp2),3-x);
                tmp3=sum(lik3,4-x); tmp3=sum(squeeze(tmp3),3-x);
                marglik1{x} = tmp1/sum(tmp1); marglik2{x} = tmp2/sum(tmp2); marglik3{x} = tmp3/sum(tmp3);
             x = 2; 
                tmp1 = sum(lik1,3-x); tmp1=sum(squeeze(tmp1),5-x);
                tmp2=sum(lik2,3-x); tmp2=sum(squeeze(tmp2),5-x);
                tmp3=sum(lik3,3-x); tmp3=sum(squeeze(tmp3),5-x);
                marglik1{x} = tmp1/sum(tmp1); marglik2{x} = tmp2/sum(tmp2); marglik3{x} = tmp3/sum(tmp3);
             x = 3; 
                tmp1 = sum(lik1,4-x); tmp1=sum(squeeze(tmp1),5-x);
                tmp2=sum(lik2,4-x); tmp2=sum(squeeze(tmp2),5-x);
                tmp3=sum(lik3,4-x); tmp3=sum(squeeze(tmp3),5-x);
                marglik1{x} = tmp1/sum(tmp1); marglik2{x} = tmp2/sum(tmp2); marglik3{x} = tmp3/sum(tmp3);                
            % plot the likelihood landscape with the maximum and true
            % value (latter for simulations only), for each individual
            % =============================================================
            
            
            %conf1
            ML(1) = myvect(p1{1}(max(marglik1{1})==marglik1{1}));
            ML(2) = myvect(p1{2}(max(marglik1{2})==marglik1{2}));
            ML(3) = myvect(p1{3}(max(marglik1{3})==marglik1{3}));
            EV(1) = sum(p1{1}(:).*marglik1{1}(:));
            EV(2) = sum(p1{2}(:).*marglik1{2}(:));
            EV(3) = sum(p1{3}(:).*marglik1{3}(:));
            [foo, poutEst{ct,1}]= RL_fitmodel_conf1(ML,data, choice1);
            [foo,  poutEst{ct,2}]= RL_fitmodel_conf1(EV,data, choice1);
            %conf2
            ML1(1) = myvect(p2{1}(max(marglik2{1})==marglik2{1}));
            ML1(2) = myvect(p2{2}(max(marglik2{2})==marglik2{2}));
            ML1(3) = myvect(p2{3}(max(marglik2{3})==marglik2{3}));
            EV1(1) = sum(p2{1}(:).*marglik2{1}(:));
            EV1(2) = sum(p2{2}(:).*marglik2{2}(:));
            EV1(3) = sum(p2{3}(:).*marglik2{3}(:));
            [foo, poutEst1{ct,1}]= RL_fitmodel_conf2(ML1,data, choice2);
            [foo,poutEst1{ct,2}]= RL_fitmodel_conf2(EV1,data, choice2);
            %conf3
            ML2(1) = myvect(p3{1}(max(marglik3{1})==marglik3{1}));
            ML2(2) = myvect(p3{2}(max(marglik3{2})==marglik3{2}));
            ML2(3) = myvect(p3{3}(max(marglik3{3})==marglik3{3}));
            EV2(1) = sum(p3{1}(:).*marglik3{1}(:));
            EV2(2) = sum(p3{2}(:).*marglik3{2}(:));
            EV2(3) = sum(p3{3}(:).*marglik3{3}(:));
            [foo, poutEst2{ct,1}]= RL_fitmodel_conf3(ML2,data, choice3);
            [foo,poutEst2{ct,2}]= RL_fitmodel_conf3(EV2,data, choice3);
            
            if plotIndividual
                % plot the grid
                figure(fhGrid);
                dims = ceil(sqrt(length(subjects)));
                subplot(dims, dims,ct);
                imagesc(lik1); hold on;
                
                [maxalpha, maxbeta] = find(lik1==max(lik1(:)));
                hline(maxalpha, {'w-', 'LineWidth', 2});
                vline(maxbeta, {'w-', 'LineWidth', 2});
                
                if simulate
                    plot(simBeta_bin*ones(1,2),[simAlpha_bin-.75 simAlpha_bin+.75],'-','color',[0.4 0.4 0.4],'linewidth',3);
                    plot([simBeta_bin-.75 simBeta_bin+.75],simAlpha_bin*ones(1,2),'-','color',[0.4 0.4 0.4],'linewidth',3);
                end
                hold off
                title(sprintf('p(data|model) all, sub %03.0f',sID));
                ytick = 0:nBin(1)/5:nBin(1);
                xtick = 0:nBin(2)/5:nBin(2);
                yticklabel = num2cell(bounds(1,1):diff(bounds(1,:)/5):bounds(1,2));
                xticklabel = num2cell(bounds(2,1):diff(bounds(2,:)/5):bounds(2,2));
                set(gca,'ytick',ytick,'yticklabel',yticklabel,'xtick',xtick,'xticklabel',xticklabel)
                ylabel(paramLabel{1}); xlabel(paramLabel{2})
                % pair 2
                figure(fhGrid1);
                dims = ceil(sqrt(length(subjects)));
                subplot(dims, dims,ct);
                imagesc(lik2); hold on;
                
                [maxalpha1, maxbeta1] = find(lik2==max(lik2(:)));
                hline(maxalpha1, {'w-', 'LineWidth', 2});
                vline(maxbeta1, {'w-', 'LineWidth', 2});
                
                if simulate
                    plot(simBeta_bin*ones(1,2),[simAlpha_bin-.75 simAlpha_bin+.75],'-','color',[0.4 0.4 0.4],'linewidth',3);
                    plot([simBeta_bin-.75 simBeta_bin+.75],simAlpha_bin*ones(1,2),'-','color',[0.4 0.4 0.4],'linewidth',3);
                end
                hold off
                title(sprintf('p(data|model) easy, sub %03.0f',sID));
                ytick = 0:nBin(1)/5:nBin(1);
                xtick = 0:nBin(2)/5:nBin(2);
                yticklabel = num2cell(bounds(1,1):diff(bounds(1,:)/5):bounds(1,2));
                xticklabel = num2cell(bounds(2,1):diff(bounds(2,:)/5):bounds(2,2));
                set(gca,'ytick',ytick,'yticklabel',yticklabel,'xtick',xtick,'xticklabel',xticklabel)
                ylabel(paramLabel{1}); xlabel(paramLabel{2})
                % pair 3
                figure(fhGrid2);
                dims = ceil(sqrt(length(subjects)));
                subplot(dims, dims,ct);
                imagesc(lik3); hold on;
                
                [maxalpha2, maxbeta2] = find(lik3==max(lik3(:)));
                hline(maxalpha2, {'w-', 'LineWidth', 2});
                vline(maxbeta2, {'w-', 'LineWidth', 2});
                

                hold off
                
                title(sprintf('p(data|model) hard, sub %03.0f',sID));
                ytick = 0:nBin(1)/5:nBin(1);
                xtick = 0:nBin(2)/5:nBin(2);
                yticklabel = num2cell(bounds(1,1):diff(bounds(1,:)/5):bounds(1,2));
                xticklabel = num2cell(bounds(2,1):diff(bounds(2,:)/5):bounds(2,2));
                set(gca,'ytick',ytick,'yticklabel',yticklabel,'xtick',xtick,'xticklabel',xticklabel)
                ylabel(paramLabel{1}); xlabel(paramLabel{2})
                
                % add the estimated choice probability to the plot
                figure(trialfh(ct)); hold on;
                [foo, foo, foo, leg] = legend;
                plot(poutEst{ct,2}.PP1s(:,1),'color',[0 .55 .55],'linewidth',2); %p(choose rew1)
                %                 leg{end+1} = 'est. p(choose rew)';
                %                 legend(leg,'location','northeastoutside');
                
                figure(trialfh1(ct)); hold on;
                [foo, foo, foo, leg] = legend;
                plot(poutEst1{ct,2}.PP2s(:,1),'color',[0 .55 .55],'linewidth',2); %p(choose rew1)
                %                 leg{end+1} = 'est. p(choose rew)';
                %                 legend(leg,'location','northeastoutside');
                
                figure(trialfh2(ct)); hold on;
                [foo, foo, foo, leg] = legend;
                plot(poutEst2{ct,2}.PP3s(:,1),'color',[0 .55 .55],'linewidth',2); %p(choose rew2)
                %                 leg{end+1} = 'est. p(choose rew)';
                %                 legend(leg,'location','northeastoutside');
                
                
            end
            
            
            % - store the marginal likelihoods across subjects, sorted by task
            % - for each parameter, compute the peak (= maximum likelihood) and
            % the mean ('expected value', i.e. sum over all bins of (value(bin)
            % * likelihood(bin))
            
            fitted.alpha.pdf(1,ct,:) = marglik1{1};
            fitted.beta.pdf(1,ct,:) = marglik1{2};
            fitted.R.pdf(1,ct,:) = marglik1{3};
            fitted.alpha1.pdf(1,ct,:) = marglik2{1};
            fitted.beta1.pdf(1,ct,:) = marglik2{2};
            fitted.R1.pdf(1,ct,:) = marglik2{3};
            fitted.alpha2.pdf(1,ct,:) = marglik3{1};
            fitted.beta2.pdf(1,ct,:) = marglik3{2};
            fitted.R2.pdf(1,ct,:) = marglik3{3};
            fitted.alpha.ml(1,ct) = ML(1);
            fitted.beta.ml(1,ct) = ML(2);
            fitted.R.ml(1,ct) = ML(3);
            fitted.alpha1.ml(1,ct) = ML1(1);
            fitted.beta1.ml(1,ct) = ML1(2);
            fitted.R1.ml(1,ct) = ML1(3);
            fitted.alpha2.ml(1,ct) = ML2(1);
            fitted.beta2.ml(1,ct) = ML2(2);
            fitted.R2.ml(1,ct) = ML2(3);
            fitted.alpha.ev(1,ct) = EV(1);
            fitted.beta.ev(1,ct) = EV(2);
            fitted.R.ev(1,ct) = EV(3);
            fitted.alpha1.ev(1,ct) = EV1(1);
            fitted.beta1.ev(1,ct) = EV1(2);
            fitted.R1.ev(1,ct) = EV1(3);
            fitted.alpha2.ev(1,ct) = EV2(1);
            fitted.beta2.ev(1,ct) = EV2(2);
            fitted.R2.ev(1,ct) = EV2(3);
            
            
        end % subject
        
        % PLOT posterior densities
        tick = 0:nBin/5:nBin;
        
        %all
        % plot posterior density of alpha across all subjects
        figure(fhParamPDF);
        subplot(2,2,1);hold on;
        plot(squeeze(fitted.alpha.pdf(:,:,:))');
        title(('alpha PDF for pair 1'))
        xticklabel = num2cell(bounds1(1,1):diff(bounds1(1,:)/5):bounds1(1,2));
        set(gca,'xtick',tick,'xticklabel',xticklabel)
        xlabel(paramLabel{1})
        ylabel('p(data|m)');
        
        % plot posterior density of beta across all subjects
        subplot(2,2,2);hold on;
        plot(squeeze(fitted.beta.pdf(:,:,:))') ;
        title(('beta PDF for pair 1'))
        xticklabel = num2cell(bounds1(2,1):diff(bounds1(2,:)/5):bounds1(2,2));
        set(gca,'xtick',tick,'xticklabel',xticklabel)
        xlabel(paramLabel{2})
        ylabel('p(data|m)');
        
        %easy
        % plot posterior density of alpha across all subjects
        figure(fhParamPDF1);
        subplot(2,2,1);hold on;
        plot(squeeze(fitted.alpha1.pdf(:,:,:))');
        title(('alpha PDF for pair 2'))
        xticklabel = num2cell(bounds2(1,1):diff(bounds2(1,:)/5):bounds2(1,2));
        set(gca,'xtick',tick,'xticklabel',xticklabel)
        xlabel(paramLabel{1})
        ylabel('p(data|m)');
        
        % plot posterior density of beta across all subjects
        subplot(2,2,2);hold on;
        plot(squeeze(fitted.beta1.pdf(:,:,:))') ;
        title(('beta PDF for pair 2'))
        xticklabel = num2cell(bounds2(2,1):diff(bounds2(2,:)/5):bounds2(2,2));
        set(gca,'xtick',tick,'xticklabel',xticklabel)
        xlabel(paramLabel{2})
        ylabel('p(data|m)');
        
        %hard
        % plot posterior density of alpha across all subjects
        figure(fhParamPDF2);
        subplot(2,2,1);hold on;
        plot(squeeze(fitted.alpha2.pdf(:,:,:))');
        title(('alpha PDF for pair 3'))
        xticklabel = num2cell(bounds3(1,1):diff(bounds3(1,:)/5):bounds3(1,2));
        set(gca,'xtick',tick,'xticklabel',xticklabel)
        xlabel(paramLabel{1})
        ylabel('p(data|m)');
        
        % plot posterior density of beta across all subjects
        subplot(2,2,2);hold on;
        plot(squeeze(fitted.beta2.pdf(:,:,:))') ;
        title(('beta PDF for pair 3'))
        xticklabel = num2cell(bounds3(2,1):diff(bounds3(2,:)/5):bounds3(2,2));
        set(gca,'xtick',tick,'xticklabel',xticklabel)
        xlabel(paramLabel{2})
        ylabel('p(data|m)');
        

        if ~plotIndividual
            
            % pair 1
            
            % plot barplots of the max and EV for alpha
            figure(fhParBar);
            
            subplot(2,2,1); barScatter(fitted.alpha.ml,[],[],true, length(subjects)); box off;
            if simulate; hold on; scatter(1:2,simpars.alpha*ones(1,2),'go','filled');hold off; end
            xlabel('all');
            ylabel('alpha (max. likelihood)');
            title('alpha - max likelihood - pair 1');
            ylim([bounds1(1,1) bounds1(1,2)])
            
            subplot(2,2,2); barScatter(fitted.alpha.ev,[],[],true,length(subjects)); box off;
            if simulate; hold on; scatter(1:2,simpars.alpha*ones(1,2),'go','filled');hold off; end
            xlabel('all');
            ylabel('alpha (expected value)');
            title('alpha - expected value - pair 1');
            ylim([bounds1(1,1) bounds1(1,2)])
            
            % plot barplots of the max and EV for beta
            subplot(2,2,3); barScatter(fitted.beta.ml,[],[],true,length(subjects)); box off;
            if simulate; hold on; scatter(1:2,simpars.beta*ones(1,2),'go','filled');hold off; end
            xlabel('all');
            ylabel('beta (max. likelihood)');
            title('beta - max likelihood - pair 1');
            ylim([bounds1(2,1) bounds1(2,2)])
            
            subplot(2,2,4); barScatter(fitted.beta.ev,[],[],true,length(subjects)); box off;
            if simulate; hold on; scatter(1:2,simpars.beta*ones(1,2),'go','filled');hold off; end
            xlabel('all');
            ylabel('beta (expected value)');
            title('beta - expected value - pair 1');
            ylim([bounds1(2,1) bounds1(2,2)])
            
            % pair 2
            
            % plot barplots of the max and EV for alpha
            figure(fhParBar1);
            
            subplot(2,2,1); barScatter(fitted.alpha1.ml,[],[],true,length(subjects)); box off;
            if simulate; hold on; scatter(1:2,simpars.alpha*ones(1,2),'go','filled');hold off; end
            xlabel('all');
            ylabel('alpha (max. likelihood)');
            title('alpha - max likelihood - pair 2');
            ylim([bounds2(1,1) bounds2(1,2)])
            
            subplot(2,2,2); barScatter(fitted.alpha1.ev,[],[],true,length(subjects)); box off;
            if simulate; hold on; scatter(1:2,simpars.alpha*ones(1,2),'go','filled');hold off; end
            xlabel('all');
            ylabel('alpha (expected value)');
            title('alpha - expected value - pair 2');
            ylim([bounds2(1,1) bounds2(1,2)])
            
            % plot barplots of the max and EV for beta
            subplot(2,2,3); barScatter(fitted.beta1.ml,[],[],true,length(subjects)); box off;
            if simulate; hold on; scatter(1:2,simpars.beta*ones(1,2),'go','filled');hold off; end
            xlabel('all');
            ylabel('beta (max. likelihood)');
            title('beta - max likelihood - pair 2');
            ylim([bounds2(2,1) bounds2(2,2)])
            
            subplot(2,2,4); barScatter(fitted.beta1.ev,[],[],true,length(subjects)); box off;
            if simulate; hold on; scatter(1:2,simpars.beta*ones(1,2),'go','filled');hold off; end
            xlabel('all');
            ylabel('beta (expected value)');
            title('beta - expected value - pair 2');
            ylim([bounds2(2,1) bounds2(2,2)])
            
            % pair 3
            % plot barplots of the max and EV for alpha
            figure(fhParBar2);
            
            subplot(2,2,1); barScatter(fitted.alpha2.ml,[],[],true,length(subjects)); box off;
            if simulate; hold on; scatter(1:2,simpars.alpha*ones(1,2),'go','filled');hold off; end
            xlabel('all');
            ylabel('alpha (max. likelihood)');
            title('alpha - max likelihood - pair 3');
            ylim([bounds3(1,1) bounds3(1,2)])
            
            subplot(2,2,2); barScatter(fitted.alpha2.ev,[],[],true,length(subjects)); box off;
            if simulate; hold on; scatter(1:2,simpars.alpha*ones(1,2),'go','filled');hold off; end
            xlabel('all');
            ylabel('alpha (expected value)');
            title('alpha - expected value - pair 3');
            ylim([bounds3(1,1) bounds3(1,2)])
            
            % plot barplots of the max and EV for beta
            subplot(2,2,3); barScatter(fitted.beta2.ml,[],[],true,length(subjects)); box off;
            if simulate; hold on; scatter(1:2,simpars.beta*ones(1,2),'go','filled');hold off; end
            xlabel('all');
            ylabel('beta (max. likelihood)');
            title('beta - max likelihood - pair 3');
            ylim([bounds3(2,1) bounds3(2,2)])
            
            subplot(2,2,4); barScatter(fitted.beta2.ev,[],[],true,length(subjects)); box off;
            if simulate; hold on; scatter(1:2,simpars.beta*ones(1,2),'go','filled');hold off; end
            xlabel('all');
            ylabel('beta (expected value)');
            title('beta - expected value - pair 3');
            ylim([bounds3(2,1) bounds3(2,2)])
            
            
            % add the trial-wise choice probabilities across subjects to the
            % probability plots
            if ~simulate
                
                % conf 1
                figure(fhData)
                
                tmp = nan(length(subjects),20);
                for isub = 1:length(subjects)
                    tmp(isub,:)= poutEst{isub,2}.PP1s(:,2); %EV only - p(rew1) for all conf
                end
                mu = mean(tmp,1);
                b = nan(20,1,1);
                b(:,1,1) = std(tmp,[],1)/sqrt(length(subjects));
                hold on
                [foo foo foo leg] = legend;
                boundedline(1:20,mu,b,'alpha','cmap',[0 .55 .55])
                leg{end+1} = 'estimated p(choice)';
                %                 legend(leg,'location','northeastoutside');
                %                 legend boxoff
            end
            
            % conf 2
            figure(fhData1)
            
            
            tmp1 = nan(length(subjects),20);
            for isub = 1:length(subjects)
                tmp1(isub,:)= poutEst1{isub,2}.PP2s(:,2); %p(rew1) for conf 1 and 2
            end
            mu1 = mean(tmp1,1);
            b1 = nan(20,1,1);
            b1(:,1,1) = std(tmp1,[],1)/sqrt(length(subjects));
            hold on
            [foo foo foo leg] = legend;
            boundedline(1:20,mu1,b1,'alpha','cmap','r')
            %             leg{end+1} = 'estimated p(choice)';
            %             legend(leg,'location','northeastoutside');
            legend boxoff
            
            
            
            % conf 3
            figure(fhData2)
            
            tmp2 = nan(length(subjects),20);
            for isub = 1:length(subjects)
                tmp2(isub,:)= poutEst2{isub,2}.PP3s(:,2); %PP (choose rew2)
            end
            mu2 = mean(tmp2,1);
            b2 = nan(20,1,1);
            b2(:,1,1) = std(tmp2,[],1)/sqrt(length(subjects));
            hold on
            [foo foo foo leg] = legend;
            boundedline(1:20,mu2,b2,'alpha','cmap','b')
            %             leg{end+1} = 'estimated p(choice)';
            %             legend(leg,'location','northeastoutside');
            legend boxoff
        end
    end
end


%==========================================================================
%% Section 4.
% Diagnostic 
%==========================================================================

ct      = 0; % reset subject counter

SS21_tot = nan (1, length(subjects));
SS21_err = nan (1, length(subjects));
r21 = nan (1, length(subjects));
SS22_tot = nan (1, length(subjects));
SS22_err = nan (1, length(subjects));
r22 = nan (1, length(subjects));
SS23_tot = nan (1, length(subjects));
SS23_err = nan (1, length(subjects));
r23 = nan (1, length(subjects));


for sID = subjects'
    for ct = 1:length(subjects)
        
       SS21_tot (ct) = nansum ((choice1_all(ct,:) - nanmean(choice1_all(ct,:))) .^2);
       SS21_err (ct) = nansum ((choice1_all(ct,:) - tmp(ct,:)) .^2);

       r21 (ct) = max (0, 1 - (SS21_err(ct) / SS21_tot(ct)));
       
       SS22_tot (ct) = nansum ((choice2_all(ct,:) - nanmean(choice2_all(ct,:))) .^2);
       SS22_err (ct) = nansum ((choice2_all(ct,:) - tmp1(ct,:)) .^2);

       r22 (ct) = max (0, 1 - (SS22_err(ct) / SS22_tot(ct)));
       
       SS23_tot (ct) = nansum ((choice3_all(ct,:) - nanmean(choice3_all(ct,:))) .^2);
       SS23_err (ct) = nansum ((choice3_all(ct,:) - tmp2(ct,:)) .^2);

       r23 (ct) = max (0, 1 - (SS23_err(ct) / SS23_tot(ct)));
    end 
end 

