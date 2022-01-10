%% VBA toolbox
%==========================================================================
% If no inputs are given, called functions will generate artificial data for both
% the learning and the test blocs and invert all data at once.
%
% IN:
%   - [data]:
%       - cues: 2 x T vector indicating the identity of the two cues
%            presented at each trial
%       - choices: 1 X T binary vector indicating if the subject chose the
%            first (choice = 0) or the second (choice = 1) cue as encoded
%            in data.cues
%       - feedbacks: 1 X T vector describing the outcome of the choice. If
%            a trial has no feedback (e.g. in test bloc), set value to NaN.
%
% OUT:
%   - posterior, out: results of the inversion
%
% /////////////////////////////////////////////////////////////////////////



subjects = [11 21 31 41 51 71 81 91 101 111 121 131 141 161 171 181  201 211 221 231 241 251 261];
rootdir = '/Users/Eliacello/Documents/PhD_DRC/RLmodel/RLtest1/';
% addpath(fullfile(rootdir,'ReinfLearn')); % modelling code
addpath(fullfile(rootdir,'util')); % some helper functions
addpath(fullfile(rootdir,'data'));% to get the feedback
addpath('/Users/Eliacello/Documents/VBA-toolbox-master/')
% code-rewrite
dataPath = fullfile(rootdir,'data'); % where our data comes from
dataPath_bv = fullfile(rootdir,'data_bv');
dataPath_ad = fullfile(rootdir,'data_ad');

% which models to run?
assymetric = 0;
simple = 0;
wsls = 0;
random_response = 0;
persevere = 0;
rew_magnitude = 1;

model_comparison = 0; 

% init
learningpos = nan (1, length(subjects));
learningneg = nan (1, length(subjects));
invtemp = nan (1, length(subjects));
chooserew = nan (1, length(subjects));
avoidpun = nan (1, length(subjects));
accmodel = nan (1, length(subjects));
baccmodel = nan (1, length(subjects));
BIC = nan (1, length(subjects));
AIC = nan (1, length(subjects));
LL = nan (1, length(subjects));
R2 = nan (1, length(subjects));
SS2_err = nan (1, length(subjects));
F = nan (1, length(subjects));
PE_all = nan (60, length(subjects));
Q_all = nan (60, length(subjects));
pred_data = nan (60, length(subjects));
%
ct      =       0; % reset subject counter

for sID = subjects
    for ct = 1:length(subjects)
        
        subTag      = sprintf('controls_S_%d_data.mat',sID(ct));% controls_S_%d_data.mat or bv_S_
        dataFile = fullfile(dataPath, subTag); %dataPath or dataPath_bv
        load(dataFile, 'data');
        alldata{ct} = data;
        
        
        data.cues = data.prep.conf';
        data.choices = data.choice';
        data.feedbacks = data.outcome';
        
        %% random:
        if random_response
            [posterior, out] = demo_random (data.choices, data.feedbacks);
        end
        
        %% win-stay-lose-shift
        if wsls
            [posterior, out] = demo_wsls_simple (data.choices, data.feedbacks);
        end
        
        %% Q-learning with two learning rates
        if assymetric
            [posterior, out, PE] = demo_QlearningAsymmetric (data);
            %             chooseA = mean(data.choices(data.prep.conf==1) == 1);
            %             avoidB = mean(data.choices (data.prep.conf==3) == 1);
            %             fprintf('Choice in test phase:\n')
            %             fprintf('  - Choose A: %3.2f%%\n',100 * chooseA);
            %             fprintf('  - Avoid  B: %3.2f%%\n',100 * avoidB);
            
            learningpos (ct) = VBA_sigmoid(posterior.muTheta(1)+posterior.muTheta(2));
            learningneg (ct) = VBA_sigmoid(posterior.muTheta(1)-posterior.muTheta(2));
            invtemp (ct) = exp(posterior.muPhi);
            
            % chooserew (ct) = chooseA;
            % avoidpun (ct) = avoidB;
            
        end
        
        %% simple Q-learning
        if simple
            [posterior, out] = demo_Qlearning (data.choices, data.feedbacks);
            %             chooseA = mean(data.choices(data.prep.conf==1) == 1);
            %             avoidB = mean(data.choices (data.prep.conf==3) == 1);
            %             fprintf('Choice in test phase:\n')
            %             fprintf('  - Choose A: %3.2f%%\n',100 * chooseA);
            %             fprintf('  - Avoid  B: %3.2f%%\n',100 * avoidB);
            
            learningpos (ct) = VBA_sigmoid(posterior.muTheta);
            invtemp (ct) = exp(posterior.muPhi);
            
        end
        
        %% Q-learning with reward magnitude
        if rew_magnitude 
            [posterior, out] = demo_QlearningR (data.choices, data.feedbacks);
            %             chooseA = mean(data.choices(data.prep.conf==1) == 1);
            %             avoidB = mean(data.choices (data.prep.conf==3) == 1);
            %             fprintf('Choice in test phase:\n')
            %             fprintf('  - Choose A: %3.2f%%\n',100 * chooseA);
            %             fprintf('  - Avoid  B: %3.2f%%\n',100 * avoidB);
            
            learningpos (ct) = VBA_sigmoid(posterior.muTheta(1));
            invtemp (ct) = exp(posterior.muPhi);
            rew_magn (ct) = VBA_sigmoid(posterior.muTheta(2));
            
            % chooserew (ct) = chooseA;
            % avoidpun (ct) = avoidB;
            
        end
        
        %% perserveration (choice kernel)
        if persevere
            [posterior, out] = demo_perseverative (data.choices, data.feedbacks);
        end
        
        
        
        %% accuracy data
        accmodel (ct) = out.fit.acc;
        baccmodel (ct) = out.fit.bacc;
        BIC (ct) = out.fit.BIC;
        AIC (ct) = out.fit.AIC;
        R2 (ct) = out.fit.R2;
        SS2_err (ct) = nansum ((data.choices - out.suffStat.gx) .^2) ;
        LL(ct) = out.fit.LL;
        F(ct) = out.F;
        %         PE_all(:,ct) = PE;
        Q_all(:, ct) = posterior.muX(1,:);
        pred_data (:,ct) = out.suffStat.gx;
        
        
        %% model comparison
        
        if model_comparison
        
        % define logEvidence_y1 = LL (model1); logEvidence_y1 = LL(model2)
        plotBayesFactor (logEvidence_y1, logEvidence_y2);
        
        % perform model selection with the VBA
        % =========================================================================
        options.verbose = false;
        
        % perform group-BMS on data generated under the full model
        [p1, o1] = VBA_groupBMC (logEvidence_y1, options);
        set (o1.options.handles.hf, 'name', 'group BMS: y_1')
        
        fprintf('Statistics in favor of the true model (m1): pxp = %04.3f (Ef = %04.3f)\n', o1.pxp(1), o1.Ef(1));
        
        % perform group-BMS on data generated under the nested model
        [p2, o2] = VBA_groupBMC (logEvidence_y2, options);
        set (o2.options.handles.hf, 'name', 'group BMS: y_2')
        
        fprintf('Statistics in favor of the true model (m2): pxp = %04.3f (Ef = %04.3f)\n', o2.pxp(2), o2.Ef(2));
        
        plotBayesFactor (logEvidence_y1, logEvidence_y2)
        [n1, x1] = VBA_empiricalDensity ((logEvidence_y1(1,:) - logEvidence_y1(2, :))');
        [n2, x2] = VBA_empiricalDensity ((logEvidence_y2(1,:) - logEvidence_y2(2 ,:))');
        hf = figure ('color' ,'w', 'name', 'demo_modelComparison: distribution of log Bayes factors');
        ha = axes ('parent', hf,'nextplot','add');
        plot (ha, x1, n1, 'color', 'r');
        plot (ha, x2, n2, 'color', 'b');
        legend (ha, {'true = model 1', 'true = model 2'});
        xlabel (ha, 'log p(y|m1) - log(y|m2)');
        ylabel (ha, 'proportion of simulations');
        
        end
    end
    
end