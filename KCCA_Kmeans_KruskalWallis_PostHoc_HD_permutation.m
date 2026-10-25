% Basic code. Including all data and no stream no subsurface analysis
% Adapted KCCA script for real pesticide data and clustering visualization
%CCA2 first line is locations 1-8. the second column is the storm number: before-0, first-1, 2- 2nd, 3- after 2
%% This version was used for the KCCA paper for the original  (not transposed) configuration.MarHeshvan 5787
% === Load your data ===
cd D:\Shulamit\PhD\Matlab_PhD\Chemistry\ML_chem
clear all 
close all
%load CCA1.txt
%load CCA1_subsurfaceTot.txt;load CCA2_Elul.txt;load CCA2_Elul_Orig.txt;x1=CCA1_subsurfaceTot;x2=CCA2_Elul;CCA2_orig=CCA2_Elul_Orig;
%load CCA1_noSubsurface.txt;load CCA2_noSubsurface.txt;load CCA2_noSubsurface_orig.txt;x1=CCA1_noSubsurface;x2=CCA2_noSubsurface;CCA2_orig=CCA2_noSubsurface_orig;cat_flow=[1 4 5 6];cat_time=[1 2 0 3];
%load CCA1_NoStream1.txt;load CCA2_NoStream1.txt;load CCA2_NoStream1_orig.txt;x1=CCA1_NoStream1;x2=CCA2_NoStream1;CCA2_orig=CCA2_NoStream1_orig;

%Best Analysis:
load CCA1_noStream.txt;load CCA2_noStream.txt;x1=CCA1_noStream;x2=CCA2_noStream;CCA2_1=CCA2_noStream;CCA2_orig=load('cca2_origNoStream.txt');
%load CCA1_noStreamWthSubsurface.txt;load CCA2_noStreamWthSubsurface.txt;load CCA2_noStreamWthSubsurface_orig.txt;x1=CCA1_noStreamWthSubsurface;x2=CCA2_noStreamWthSubsurface;CCA2_orig=CCA2_noStreamWthSubsurface_orig;
%load CCA1_subsurfaceTot.txt;load CCA2_subsurfaceTot.txt;load CCA2_subsurfaceTot_orig.txt;x1=CCA1_subsurfaceTot;x2=CCA2_subsurfaceTot;CCA2_orig=CCA2_subsurfaceTot_orig;%In this set manhole and pipe outlet where unified as one flowpaths since1. it is the same flowpats. 2. KW is sensitive to group size

%NoSoil
%load CCA1_noSoil.txt;load CCA2_noSoil.txt;x1=CCA1_noSoil;x2=CCA2_noSoil;CCA2_1=CCA2_noSoil;CCA2_orig=load('CCA2_orig_NoSoil.txt');


%load GW1.txt;load GW2.txt;x1=GW1;x2=GW2;GW2_1=GW2;
%for runoff
%CCA2=CCA2(:,1);In a case we want to isolate a controlling factor
%load CCA1_shrt.txt;load CCA2.txt;x1=CCA1_shrt;x2=CCA2;CCA2_1=CCA2;
%CCA2=CCA2(:,2);x1=CCA1_shrt;x2=CCA2;%In a case we want to isolate a controlling factor
% load CCA1.txt
%load CCA1_mob.txt;load CCA2_mob.txt;x1=CCA1_mob;x2=CCA2_mob;CCA2_1=CCA2_mob;%For mobility parameters the matrix was transposed, since features in CCA2 are on the pesticides and not sampling locations
%load CCA1_mob_rnff.txt;load CCA2_mob_rnff.txt;x1=CCA1_mob_rnff;I=find(CCA1_mob_rnff==0);CCA1(I)=0.0025;x2=CCA2_mob_rnff;CCA2_1=CCA2_mob_rnff;%For mobility parameters the matrix was transposed, since features in CCA2 are on the pesticides and not sampling locations

%load CCA1_shrt_mob.txt;load CCA2_mob.txt;x1=CCA1_shrt_mob;x2=CCA2_mob;CCA2_1=CCA2_mob;
%CCA2_mob=CCA2_mob(:,4);x1=CCA1_shrt_mob;x2=CCA2_mob;
%load CCA1_shrt.txt;load CCA2.txt;
%%
%Creating scaling of data not to have range.
% epsilon = 1e-4;

% Log-transform the data
% X1_log = log(CCA1_shrt + epsilon);
% X2_log = log(CCA2 + epsilon);
% 
% % Normalize (z-score)
% X1_norm = (X1_log - mean(X1_log)) ./ std(X1_log);
% X2_norm = (X2_log - mean(X2_log)) ./ std(X2_log);
% x1 = X1_norm; % Assign normalized data to x1
% x2 = x2_norm; % Assign normalized data to x2


%%
%A=CCA1_shrt(:,1);
% shuffled_A = A(randperm(length(A)));
% CCA1_shrt(:,1)=shuffled_A;
%filename='loc_noStream.txt';
%filename='loc_noSubsurface.txt';
%filename='loc_name.txt';

%Best analysis
filename='loc_nameNoStream.txt';
%filename='loc_nameNoSoil.txt';
%filename='loc_noStreamWthSubsurface'%Also for tot

%filename='nameShrt.txt';
%filename='GW_name.txt';
T=readtable(filename);
name=T{:,1};
% Replace these with your actual data matrices
%x1 = CCA1_shrt;%pesticide_data_matrix;  % samples x pesticides
% x1=CCA1;
% x2 = CCA2;%environmental_matrix;   % samples x features or metadata

labels = name; % vector of sample types (numeric or categorical)

% === Preprocessing ===
N = size(x1,1);

% Centering
x1 = x1 - mean(x1);
% x2 = x2 - mean(x2);x2 represents categorical controls (flowpath, timing, field) — these must remain integers or one-hot encoded values, never normalized.

% Normalize variance
x1 = x1 * sqrt(diag(1./diag(x1'*x1)));
% x2 = x2 * sqrt(diag(1./diag(x2'*x2)));x2 represents categorical controls (flowpath, timing, field) — these must remain integers or one-hot encoded values, never normalized.

% === Parameters ===
kerneltype = 'gauss';
%kernelpar = 1; %sigma     % Tune this later
D = pdist2(x1, x1);
kernelpar = median(D(:)); % Good default
numeig=3
reg = 1e-2;         % Regularization
%X1=x1;X2=x2;
% === Run KCCA ===
[y1, y2, beta, G1, G2,alpha1,alpha2] = km_kcca(x1, x2, kerneltype, kernelpar, reg,numeig);
beta
for k = 1:3
    fprintf('Component %d: beta = %.4f, corr = %.4f\n', ...
        k, beta(k), corr(y1(:,k),y2(:,k)));
end
% === Clustering and Visualization ===
Z = [y1, y2];  % 2D embedding

% Plot labeled groups
figure;
gscatter(Z(:,1), Z(:,4), labels);% Z(:,[1 4]) because numeign=3 means three first are by x1 and three last, i.e. 4:6, are by x2. Best correlated are the first ones in each set (i.e.columns 1 and 4).
% (Z(:,1), Z(:,2), labels);Scatter is based on the vectors in the kernel space that are most correlated. 
xlabel('KCCA 1'); ylabel('KCCA 2');
title('KCCA Projection by Sample Type');

R=corr(x1,y1);%This function checks how the x1 original data are correlated with the latent space of y1.
%i.e. how they most affect the y1 vector. Indeed check which pesticides
%most effect the y1 in the latent space, thus the KCCA
% Optional: run unsupervised clustering (k-means)

 
%% 
%For coloring according to location:
% Assume y1 and y2 are 107×1 (or more) vectors from KCCA
figure (2)
siteCategory = CCA2_orig(:,1);%CCA2_orig(:,1); % Extract site categories from CCA2
% siteCategory=[8	8	2	3	5	5	6	6	6	7	7	5	5	5 ...
%               5	6	8	3	3	5	5	5	6	6	7	7	7	7 ...
%           	5	5   6	6	5	5	6	6	7	8	1	1	1	1 ...
%         	1	1	1	1   1	1	1	1	1	1	1	1	1	1 ...
%         	1	1	1	1	1	1   8	8	8	8	8	8	8	8 ...
%         	8	8	8	8	8	2	2	2   2	2	2	8	8	8 ...
%         	8	8	8	4	4	4	4	4	4	4   4	4	4	4 ...
%         	4	4	4	4	4	4	4	4	4];
% siteCategory=[2	3	5	5	6	6	6	7	7	5	5	5 ...
%               5	6	3	3	5	5	5	6	6	7	7	7 ...
%           	  7 5	5   6	6	5	5	6	6	7	1	1 ...
%           	  1 1   1	1	1	1   1	1	1	1	1	1 ...	
%               1	1	1	1   1	1	1	1	1	1   2	2 ...
%           	  2 2	2	2	4	4	4	4	4	4	4   4 ...
%           	  4	4	4   4	4	4	4	4	4	4	4	4];%Only no
%           	  stream
% siteCategory=[3	3	4	4	4	5	5	3	3	3   3	4	3	3	3	4	4	5	5	5  5 3	3   4	4	3	3	4	4	5   1   1   1   1   1	1	1 ...	
%               1 1	1   1   1	1	1   1	1	1	1   1	1   1   1   1 ...
%               1   2	2	2  2 2	2   2   2   2	2	2   2	2	2	2 ...
%               2	2	2	2	2];%NoStream, No subsurface: manhole+pipe (compare 'name')
% ;%CCA2_1(:,1);
gscatter(y1(:,1), y2(:,1), siteCategory)
%legend('Location', 'best')  % will now show only 8 groups
legend({'GW','S','OF','SD','PD'}, 'Location', 'best')
%legend({'GW','PI','M','S','OF','SC','PC'},'Location','best')%GW-groundwater;PI-pipe;M-manhole; S-soil;OF-overland; flow;SC-secondary channel;promary channel
%legend({'GW','OF','SC','PC'}, 'Location', 'best')%No Soil

xlabel('KCCA projection 1')
ylabel('KCCA projection 2')
title('KCCA Results Colored by Site Category')
%If you want to customize colors or use scatter for finer control:

% colors = lines(8); % 8 distinguishable colors
% hold on
% for i = 1:8
%     idx = siteCategory == i;
%     scatter(y1(idx,1), y2(idx,1), 60, colors(i,:), 'filled', 'DisplayName', ['Group ' num2str(i)])
% end
% legend('Location', 'best')
% hold off
%% Kmeans
 k = 5;%numel(unique(labels));8
 [idx, ~] = kmeans(Z, k);
% figZ(:,figure('name','Kmeans');
% gscatter(1), Z(:,2), idx);
% xlabel('KCCA 1'); ylabel('KCCA 2');
% title('KCCA Clustering (k-means)');

% === Report ===
%fprintf('Canonical correlation (first component): %f\n', beta);


%% Kruskal–Wallis + violin plot See chat:"correlation impact on concentration"
% =========================================
%1. Load data
%=========================================
cca2_orig =CCA2_orig;%CCA2_noSubsurface_orig;% load('cca2_origNoStream.txt');   % columns: [flowpath, timing, location]
% Y1 must already be in workspace (107 x m), canonical variates from KCCA

%% ===================== CONFIG & INPUTS =====================
% Assumed available:
%   y1 : [n x k] canonical variates (use y1(:,1))
%   x2 : [n x 3] integer labels:
%        x2(:,1) = Flowpath (1..7)
%        x2(:,2) = Timing   (1..4)
%        x2(:,3) = Field    (1..2)  -> 1=East, 2=West  (edit if opposite)

scores = y1(:,1);                 % canonical variate to analyze

flow_id  = cca2_orig(:,1);               % 1..7
time_id  = cca2_orig(:,2);               % 1..4
field_id = cca2_orig(:,3);               % 1..2 (East/West)

% Names (exactly as requested)
%flow_names  = cellstr(["Groundwater","Overland flow","Secondary trench", "Primary trench"]);d=4;% convert to cell array for strjoin
%flow_names  = cellstr(["Groundwater","Pipe","Manhole","Overland flow","Primary trench","Secondary trench"]);d=6;cat_flow=[1 2 3 4 5 6];cat_time=[1 2 3 4];% convert to cell array for strjoin
%flow_names  = cellstr(["Groundwater","Pipe","Manhole","Soil","Overland flow","Primary trench","Secondary trench", "Stream"]);d=8;% convert to cell array for strjoin

%Best analysis
%flow_names  = cellstr(["Groundwater","Soil","Overland flow","Secondary trench","Primary trench"]);d=5;cat_flow=[1 2 3 4 5];cat_time=[1 2 3 4];%No_stream=no stream+subsurface, both subsurface at the other edge of the linear
flow_names  = cellstr(["GW","S","OF","SD","PD"]);d=5;cat_flow=[1 2 3 4 5];cat_time=[1 2 3 4];%No_stream=no stream+subsurface, both subsurface at the other edge of the linear
%flow_names  = cellstr(["GW","OF","SC","PC"]);d=4;cat_flow=[1 2 3 4 5];cat_time=[1 2 3 4];%No_stream=no stream+subsurface, both subsurface at the other edge of the linear. No SOil
%flow_names  = cellstr(["Groundwater","Pipe","Manhole","Soil","Overland flow","Secondary trench","Primary trench"]);d=7;cat_flow=1:7;cat_time=[1 2 3 4];%Stream with subNo_stream=no stream+subsurface, both subsurface at the other edge of the linear
%flow_names  = cellstr(["GW","P","M","S","OF","SD","PD"]);d=7;cat_flow=1:7;cat_time=[1 2 3 4];%Stream with subNo_stream=no stream+subsurface, both subsurface at the other edge of the linear
%flow_names  = cellstr(["Groundwater","Subsurface","Soil","Overland flow","Secondary trench","Primary trench"]);d=6;cat_flow=1:6;cat_time=[1 2 3 4];%%Stram with total sub surface. No_stream=no stream+subsurface, both subsurface at the other edge of the linear
time_names  = cellstr(["1st storm","2nd storm","Before 2nd storm","After 2nd storm"]);
%time_names  = cellstr(["1st storm","2nd storm","Before 2nd storm","After 2nd storm"]);
field_names = cellstr(["East","West"]);     % adjust order if needed
 
% Categorical groupings (keeps numeric order on the axis)
gFlow  = categorical(flow_id , cat_flow, flow_names);%If including manhole+pipes+stream: 1:8
gTime  = categorical(time_id , cat_time, time_names);
gField = categorical(field_id, 1:2, field_names);

%% ===================== PLOTTING =====================
figure('Units','normalized','Position',[0.05 0.15 0.9 0.7]);

subplot(1,3,1);
plot_violin_with_counts(scores, gFlow, 'Distribution per Flowpath', 'Flowpath');
ylabel('Canonical variate y1(:,1)');

subplot(1,3,2);
plot_violin_with_counts(scores, gTime, 'Distribution per Timing', 'Timing');

subplot(1,3,3);
plot_violin_with_counts(scores, gField, 'Distribution per Field', 'Field');

%% ===================== STATS: KW + POST-HOC =====================
fprintf('\n==== Kruskal–Wallis & Post-hoc (Flowpath) ====\n');
[p_kw, tbl, stats]=kw_posthoc(scores, flow_id, flow_names);

fprintf('\n==== Kruskal–Wallis & Post-hoc (Timing) ====\n');
[p_kw, tbl, stats]=kw_posthoc(scores, time_id, time_names);

fprintf('\n==== Kruskal–Wallis & Post-hoc (Field) ====\n');
[p_kw, tbl, stats]=kw_posthoc(scores, field_id, field_names);
%% === Hellinger Distance analysis (HD) + permutation test ===
fprintf('\n==== Hellinger Distance (HD) — Flowpath groups ====\n');

groups = flow_id;
uniqueGroups = unique(groups);
nG = numel(uniqueGroups);
HD = nan(nG);
HD(1:nG+1:end) = 0;

% IMPORTANT:
% The same bin edges are used for the observed and all permuted datasets.
% The permutation test shuffles FLOWPATH LABELS only, thereby preserving:
%   1) the y1 canonical scores,
%   2) the exact group sizes,
%   3) the overall score distribution.
% It tests whether the observed HD is larger than expected from random
% assignment of these scores to flowpath groups.
%
% This is a post-hoc permutation test on the fixed KCCA scores. It does NOT
% refit KCCA in every permutation.

nBins = 29;   % 30 edges -> 29 bins, as in the previous version
bins = linspace(min(scores), max(scores), nBins+1);

%% --------------------- Observed pairwise HD ---------------------
for i = 1:nG
    for j = i+1:nG
        Pi = histcounts(scores(groups==uniqueGroups(i)), bins, ...
            'Normalization','probability');
        Pj = histcounts(scores(groups==uniqueGroups(j)), bins, ...
            'Normalization','probability');

        HD(i,j) = (1/sqrt(2)) * sqrt(sum((sqrt(Pi)-sqrt(Pj)).^2));
        HD(j,i) = HD(i,j);
    end
end

% visualize observed HD
figure('Name','Hellinger Distance (Flowpath)', ...
    'Units','normalized','Position',[0.2 0.2 0.35 0.5]);
imagesc(HD,[0 1]); colorbar;
title('Observed Hellinger distance between flowpath distributions');
xticks(1:nG); yticks(1:nG);
xticklabels(flow_names(1:nG)); yticklabels(flow_names(1:nG)); xtickangle(30);

disp('Observed pairwise HD matrix:');
disp(array2table(HD, ...
    'VariableNames',flow_names(1:nG), ...
    'RowNames',flow_names(1:nG)));

%% --------------------- Permutation test for HD ---------------------
nPerm = 5000;
rng(1,'twister');   % reproducible permutations

% Store only upper-triangle pairwise distances.
pairMask = triu(true(nG),1);
[pairI,pairJ] = find(pairMask);
nPairs = numel(pairI);

HD_obs_pairs  = zeros(nPairs,1);
HD_perm_pairs = zeros(nPerm,nPairs);

for q = 1:nPairs
    HD_obs_pairs(q) = HD(pairI(q),pairJ(q));
end

fprintf('Running %d label permutations for Hellinger distance...\n',nPerm);

for b = 1:nPerm
    permGroups = groups(randperm(numel(groups)));  % preserves group sizes

    for q = 1:nPairs
        i = pairI(q);
        j = pairJ(q);

        Pi = histcounts(scores(permGroups==uniqueGroups(i)), bins, ...
            'Normalization','probability');
        Pj = histcounts(scores(permGroups==uniqueGroups(j)), bins, ...
            'Normalization','probability');

        HD_perm_pairs(b,q) = (1/sqrt(2)) * ...
            sqrt(sum((sqrt(Pi)-sqrt(Pj)).^2));
    end
end

% Empirical one-sided p-value:
% probability under random labels of HD >= observed HD.
pHD_pairs = zeros(nPairs,1);
for q = 1:nPairs
    pHD_pairs(q) = (1 + sum(HD_perm_pairs(:,q) >= HD_obs_pairs(q))) ...
                    / (nPerm + 1);
end

% Family-wise error control across all flowpath pairs using a max-statistic.
% For every permutation, take the largest HD among all pairwise comparisons.
maxHD_perm = max(HD_perm_pairs,[],2);
pHD_FWER_pairs = zeros(nPairs,1);
for q = 1:nPairs
    pHD_FWER_pairs(q) = (1 + sum(maxHD_perm >= HD_obs_pairs(q))) ...
                         / (nPerm + 1);
end

% Null distribution summaries
HD_null_mean = mean(HD_perm_pairs,1)';
HD_null_sd   = std(HD_perm_pairs,0,1)';
HD_null_q95  = prctile(HD_perm_pairs,95,1)';

% Put results back into symmetric matrices
pHD      = nan(nG);
pHD_FWER = nan(nG);
HDnull95 = nan(nG);

for q = 1:nPairs
    i = pairI(q);
    j = pairJ(q);

    pHD(i,j) = pHD_pairs(q);
    pHD(j,i) = pHD_pairs(q);

    pHD_FWER(i,j) = pHD_FWER_pairs(q);
    pHD_FWER(j,i) = pHD_FWER_pairs(q);

    HDnull95(i,j) = HD_null_q95(q);
    HDnull95(j,i) = HD_null_q95(q);
end

% Pairwise result table
Group1 = strings(nPairs,1);
Group2 = strings(nPairs,1);
for q = 1:nPairs
    Group1(q) = string(flow_names{pairI(q)});
    Group2(q) = string(flow_names{pairJ(q)});
end

T_HDperm = table(Group1,Group2,HD_obs_pairs,HD_null_mean,HD_null_sd, ...
    HD_null_q95,pHD_pairs,pHD_FWER_pairs, ...
    'VariableNames',{'Group1','Group2','ObservedHD','NullMeanHD', ...
    'NullSD','Null95pct','PermutationP','MaxStatFWER_P'});

disp(' ');
disp('==== Pairwise Hellinger permutation results ====');
disp(T_HDperm);

%% --------------------- Permutation p-value maps ---------------------
% Raw empirical permutation p-values
figure('Name','HD permutation p-values (Flowpath)', ...
    'Units','normalized','Position',[0.2 0.2 0.35 0.5]);
imagesc(-log10(pHD)); colorbar;
title('-log_{10}(p): HD label-permutation test');
xticks(1:nG); yticks(1:nG);
xticklabels(flow_names(1:nG)); yticklabels(flow_names(1:nG)); xtickangle(30);

% Max-statistic FWER-corrected p-values
figure('Name','HD permutation FWER p-values (Flowpath)', ...
    'Units','normalized','Position',[0.2 0.2 0.35 0.5]);
imagesc(-log10(pHD_FWER)); colorbar;
title('-log_{10}(p): HD permutation, max-stat FWER');
xticks(1:nG); yticks(1:nG);
xticklabels(flow_names(1:nG)); yticklabels(flow_names(1:nG)); xtickangle(30);

%% --------------------- Global HD permutation test ---------------------
% One global statistic: mean of all pairwise HD values.
HD_global_obs  = mean(HD_obs_pairs);
HD_global_perm = mean(HD_perm_pairs,2);
pHD_global = (1 + sum(HD_global_perm >= HD_global_obs)) / (nPerm + 1);

fprintf('\nGlobal mean pairwise HD = %.4f\n',HD_global_obs);
fprintf('Permutation null mean   = %.4f\n',mean(HD_global_perm));
fprintf('Global HD permutation p = %.6g\n',pHD_global);
fprintf(['Interpretation: pairwise HD is descriptive effect size; ', ...
         'PermutationP tests each pair against random flowpath labels; ', ...
         'MaxStatFWER_P controls family-wise error across all flowpath pairs.\n']);


meanY1 = grpstats(scores, flow_id, 'mean');

% This gives: average Y₁ per flowpath

%% ===================== HELPERS =====================
function plot_violin_with_counts(vals, groups, ttl, xlab)
    % One call builds all violins with labels + n per group
    violinplot(vals, groups);
    title(ttl); xlabel(xlab); grid on; xtickangle(20);

    % annotate sample sizes above violins
    hold on
    [grpIdx, grpCats] = findgroups(groups);
    nPer = splitapply(@numel, vals, grpIdx);
    yl = ylim; ytxt = yl(2) - 0.04*range(yl);
    for i = 1:numel(nPer)
        text(i, ytxt, sprintf('n=%d', nPer(i)), 'HorizontalAlignment','center');
    end
    hold off

    % Print the x-axis order to console
    disp([ttl ' — axis order:']); disp(string(categories(groups)).');
end

function [p_kw, tbl, stats]=kw_posthoc(vals, ids, names)
    % Global Kruskal–Wallis
    [p_kw, tbl, stats] = kruskalwallis(vals, ids, 'off');
    H = tbl{2,4};   % KW statistic
    fprintf('Global KW p-value: %.3g\n', p_kw);
    fprintf('Global KW H: %.3g\n', H);
    % fprintf('Global KW stats: %.3g\n', stats);
    % MATLAB performs these steps internally:
% Convert all values into ranks across the full dataset.
% Compute H statistic using:
% Adjust for ties (using a correction factor if tied ranks exist).
% Compute the p-value by comparing H to the chi-square distribution with 
% where chi2cdf is the cumulative χ² distribution function.
% So your MATLAB output:
% Global KW p-value: 1.16e-13
% is just the p = 1 - χ²cdf(H, k−1) result — MATLAB hides H unless you inspect the second output (tbl).

    % Post-hoc pairwise (Dunn–Šidák controls FWER; Bonferroni also OK)
    res = multcompare(stats, 'CType','dunn-sidak', 'Display','off');
    % res columns: [grp1 grp2 lowerCI diff upperCI pValue]

    T = table( ...
        names(res(:,1))', names(res(:,2))', res(:,4), res(:,6), ...
        'VariableNames', {'Group1','Group2','Difference','pValue'});
    disp(T);

    % Optional: significance matrix (-log10 p) for a quick glance
    u = unique(ids); k = numel(u);
    P = nan(k); 
    for r = 1:size(res,1)
        i = res(r,1); j = res(r,2);
        P(i,j) = res(r,6); P(j,i) = res(r,6);
    end
    figure('Name',['Post-hoc map — ' strjoin(names,' / ')],'Units','normalized','Position',[0.2 0.2 0.35 0.5]);
    imagesc(-log10(P)); colorbar; 
    title('-log_{10}(p) pairwise'); 
    xticks(1:k); yticks(1:k);
    xticklabels(names); yticklabels(names); xtickangle(30);
end
%% Debugging
disp('--- Debug info ---');
disp(['size(scores) = ', mat2str(size(scores))]);
disp(['size(flow_id) = ', mat2str(size(flow_id))]);
disp(['class(flow_id) = ', class(flow_id)]);
disp('first 10 values of flow_id:');
disp(flow_id(1:10));

