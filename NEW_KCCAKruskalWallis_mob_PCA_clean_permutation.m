%% This version was used for the KCCA paper for the transposed configuration.MarHeshvan 5787
%  Adapted KCCA script for real pesticide data and clustering visualizationKoc_
%CCA2 first line is locations 1-8. the second column is the storm number: before-0, first-1, 2- 2nd, 3- after 2
%NEW...._clean and NEW...._PCA are exactly the same (accept the 'idxKoc' definition).
% === Load your data ===
cd D:\Shulamit\PhD\Matlab_PhD\Chemistry\ML_chem
clear all 
close all
clc
%load CCA1.txt
%load CCA1.txt;load CCA2.txt;x1=CCA1;x2=CCA2;CCA2_1=CCA2;
%load CCA1_noStream.txt;load CCA2_noStream.txt;x1=CCA1_noStream;x2=CCA2_noStream;CCA2_1=CCA2_noStream;
%load GW1.txt;load GW2.txt;x1=GW1;x2=GW2;GW2_1=GW2;
%for runoff
%CCA2=CCA2(:,1);In a case we want to isolate a controlling factor
%load CCA1_shrt.txt;load CCA2.txt;x1=CCA1_shrt;x2=CCA2;CCA2_1=CCA2;
%CCA2=CCA2(:,2);x1=CCA1_shrt;x2=CCA2;%In a case we want to isolate a controlling factor
% load CCA1.txt

%MOB (ones with %% give high, not significant p values)
%%load CCA1_mob.txt;load CCA25_mob.txt;load CCA25_mob_orig.txt;x1=CCA1_mob;x2=CCA25_mob;CCA2_1=CCA25_mob;Orig_CCA2_mob=CCA25_mob_orig;%For mobility parameters the matrix was transposed, since features in CCA2 are on the pesticides and not sampling locations
%%load CCA1_pca_mob.txt;load CCA25_mob.txt;load CCA25_mob_orig.txt;x1=CCA1_pca_mob;x2=CCA25_mob;Orig_CCA2_mob=CCA25_mob_orig;CCA2_mob=CCA25_mob;CCA1_mob=CCA1_pca_mob;%removing flowpath that are less significant regarding the data distribution. cca2 is same as full data (see upper line). Including all compounds that we could find the koc, Koc, solubility half-life and dominant flowpats
%load CCA11_mob.txt;load CCA25_mob.txt;load CCA25_mob_orig.txt;x1=CCA11_mob;x2=CCA25_mob;CCA2_1=CCA25_mob;Orig_CCA2_mob=CCA25_mob_orig;%For mobility parameters the matrix was transposed, since features in CCA2 are on the pesticides and not sampling locations
%load CCA1_full_mob.txt;load CCA25_mob.txt;load CCA25_mob_orig.txt;Orig_CCA2_mob=CCA25_mob_orig;x1=CCA1_full_mob;x2=CCA25_mob;CCA2_mob=CCA25_mob;CCA1_mob=CCA1_full_mob;orig_CCA2=CCA25_mob_orig;%including all compounds that we could find the koc, Koc, solubility half-life and dominant flowpats
%%load CCA1_full_full_mob.txt;load CCA25_mob.txt;load CCA25_mob_orig.txt;Orig_CCA2_mob=CCA25_mob_orig;x1=CCA1_full_full_mob;x2=CCA25_mob;CCA2_mob=CCA25_mob;CCA1_mob=CCA1_full_full_mob;orig_CCA2=CCA25_mob_orig;%including all compounds that we could find the koc, Koc, solubility half-life and dominant flowpats, including all p2-9 and all GW
%load CCA1_pca_mob.txt;load CCA2_full_mob.txt;load CCA25_mob_orig.txt;x1=CCA1_pca_mob;x2=CCA2_full_mob;CCA2_mob=CCA2_full_mob;CCA1_mob=CCA1_pca_mob;Orig_CCA2_mob=CCA25_mob_orig;%removing flowpath that are less significant regarding the data distribution. cca2 is same as full data (see upper line). Including all compounds that we could find the koc, Koc, solubility half-life and dominant flowpats
%%load CCA1_pca_mob.txt;load CCA25_mob.txt;load CCA25_mob_orig.txt;x1=CCA1_pca_mob;x2=CCA25_mob;CCA2_mob=CCA25_mob;Orig_CCA2_mob=CCA25_mob_orig;CCA1_mob=CCA1_pca_mob;%removing flowpath that are less significant regarding the data distribution. cca2 is same as full data (see upper line). Including all compounds that we could find the koc, Koc, solubility half-life and dominant flowpats

%load CCA1_shrt_mob.txt;load CCA2_shrt_mob.txt;x1=CCA1_shrt_mob;x2=CCA2_shrt_mob;CCA2_mob=CCA2_shrt_mob;CCA1_mob=CCA1_shrt_mob;%including full data set (almost) compounds that we could find the koc, Koc, solubility half-life and dominant flowpats
%load CCA1_mob_flwpth.txt;load CCA25_mob.txt;load CCA25_mob_orig.txt;x1=CCA1_mob_flwpth;x2=CCA25_mob_orig;CCA2_1=CCA25_mob_orig;Orig_CCA2_mob=CCA25_mob_orig;CCA1_mob=CCA1_mob_flwpth;%Including east GW no PW14For mobility parameters the matrix was transposed, since features in CCA2 are on the pesticides and not sampling locations
%%load CCA11_mob.txt;load CCA25_mob.txt;load CCA25_mob_orig.txt;x1=CCA11_mob;x2=CCA25_mob;CCA2_1=CCA25_mob;Orig_CCA2_mob=CCA25_mob_orig;CCA1_mob=CCA11_mob;%Including PW14 and east GW
%load CCA1_selected.txt;load CCA25_mob.txt;load CCA25_mob_orig.txt;x1=CCA1_selected;x2=CCA25_mob;CCA2_1=CCA25_mob;Orig_CCA2_mob=CCA25_mob_orig;CCA1_mob=CCA1_selected;%Including PW14 and east GW
%%load CCA12_mob.txt;load CCA25_mob.txt;load CCA25_mob_orig.txt;x1=CCA12_mob;x2=CCA25_mob;CCA2_1=CCA25_mob;Orig_CCA2_mob=CCA25_mob_orig;CCA1_mob=CCA12_mob;%Including PW14 no east GW
%%load CCA13_mob.txt;load CCA25_mob.txt;load CCA25_mob_orig.txt;x1=CCA13_mob;x2=CCA25_mob;CCA2_1=CCA25_mob;Orig_CCA2_mob=CCA25_mob_orig;CCA1_mob=CCA13_mob;%Including PW14 no east GW
%%load CCA14_mob.txt;load CCA25_mob.txt;load CCA25_mob_orig.txt;x1=CCA14_mob;x2=CCA25_mob;CCA2_1=CCA25_mob;Orig_CCA2_mob=CCA25_mob_orig;CCA1_mob=CCA14_mob;%Including west field only, no PW14 
%%load CCA15_mob.txt;load CCA25_mob.txt;load CCA25_mob_orig.txt;x1=CCA15_mob;x2=CCA25_mob;CCA2_1=CCA25_mob;Orig_CCA2_mob=CCA25_mob_orig;CCA1_mob=CCA15_mob;%Including west field only, with PW14
%load CCA1_no_subsurface_mob.txt;load CCA25_mob.txt;load CCA25_mob_orig.txt;x1=CCA1_no_subsurface_mob;x2=CCA25_mob;CCA2_1=CCA25_mob;Orig_CCA2_mob=CCA25_mob_orig;CCA1_mob=CCA1_no_subsurface_mob;%Including west field only, with PW14
%%load CCA1_no_stream_mob.txt;load CCA25_mob.txt;load CCA25_mob_orig.txt;x1=CCA1_no_stream_mob;x2=CCA25_mob;CCA2_1=CCA25_mob;Orig_CCA2_mob=CCA25_mob_orig;CCA1_mob=CCA1_no_stream_mob;%Inc

%load CCA1_mob_noStrmNoSub.txt;load CCA2_mob_noStrmNoSub.txt;load CCA2_mob_noStrmNoSub_orig.txt;x1=CCA1_mob_noStrmNoSub;x2=CCA2_mob_noStrmNoSub;CCA2_1=CCA2_mob_noStrmNoSub;CCA25_mob_orig=CCA2_mob_noStrmNoSub_orig;Orig_CCA2_mob=CCA2_mob_noStrmNoSub_orig;CCA1_mob=CCA1_mob_noStrmNoSub;%Inc
%load CCA1_mob_noStrmNoSub.txt;load CCA2_mob_noStrmNoSub1.txt;load CCA2_mob_noStrmNoSub_orig1.txt;x1=CCA1_mob_noStrmNoSub;x2=CCA2_mob_noStrmNoSub1;CCA2_1=CCA2_mob_noStrmNoSub1;CCA25_mob_orig=CCA2_mob_noStrmNoSub_orig1;Orig_CCA2_mob=CCA2_mob_noStrmNoSub_orig1;CCA1_mob=CCA1_mob_noStrmNoSub;%Inc
%load CCA1_Red_mob_noStrmNoSub11.txt;load CCA2_mob_noStrmNoSub11.txt;load CCA2_mob_noStrmNoSub_orig11.txt;x1=CCA1_Red_mob_noStrmNoSub11;x2=CCA2_mob_noStrmNoSub11;CCA2_1=CCA2_mob_noStrmNoSub11;CCA25_mob_orig=CCA2_mob_noStrmNoSub_orig11;Orig_CCA2_mob=CCA2_mob_noStrmNoSub_orig11;CCA1_mob=CCA1_Red_mob_noStrmNoSub11;
%Short compounds series
%%load CCA13_mob.txt;load CCA25_mob.txt;load CCA25_mob_orig.txt;x1=CCA13_mob;x2=CCA25_mob;CCA2_1=CCA25_mob;Orig_CCA2_mob=CCA25_mob_orig;CCA1_mob=CCA13_mob;
%%load CCA16_mob.txt;load CCA2_mob_shrt.txt;load CCA2_mob_shrt_orig.txt;x1=CCA16_mob;x2=CCA2_mob_shrt;CCA2_1=CCA2_mob_shrt;Orig_CCA2_mob=CCA2_mob_shrt_orig;CCA1_mob=CCA16_mob;%average concentration per flowpaths
%%load CCA17_mob.txt;load CCA2_mob_shrt.txt;load CCA2_mob_shrt_orig.txt;x1=CCA17_mob;x2=CCA2_mob_shrt;CCA2_1=CCA17_mob;Orig_CCA2_mob=CCA2_mob_shrt_orig;CCA1_mob=CCA17_mob;%reduced matrix for stability. Only selected by PCA and some more manually were rmoved

%Best results
%load CCA1_mob_noStrmNoSub.txt;load CCA2_mob_noStrmNoSub.txt;load CCA2_mob_noStrmNoSub.txt;load CCA2_mob_noStrmNoSub_orig.txt;x1=CCA1_mob_noStrmNoSub;x2=CCA2_mob_noStrmNoSub;CCA2_1=CCA2_mob_noStrmNoSub;CCA25_mob_orig=CCA2_mob_noStrmNoSub_orig;Orig_CCA2_mob=CCA2_mob_noStrmNoSub_orig;CCA1_mob=CCA1_mob_noStrmNoSub;%Inc

%Elul
%load CCA1_mob_noStrmNoSub.txt;load CCA2_mob_noStrmNoSub_Elul.txt;load CCA2_mob_noStrmNoSub_Elul.txt;load CCA2_mob_noStrmNoSub_orig.txt;x1=CCA1_mob_noStrmNoSub;x2=CCA2_mob_noStrmNoSub_Elul;CCA2_1=CCA2_mob_noStrmNoSub_Elul;%CCA25_mob_orig=CCA2_mob_noStrmNoSub_orig;Orig_CCA2_mob=CCA2_mob_noStrmNoSub_orig;CCA1_mob=CCA1_mob_noStrmNoSub;%Inc
%load CCA1_no_sub_mob.txt;load CCA2_no_sub_mob.txt;load CCA2_no_sub_orig_mob.txt;x1=CCA1_no_sub_mob;x2=CCA2_no_sub_mob;CCA2_1=CCA1_no_sub_mob;CCA25_mob_orig=CCA2_no_sub_orig_mob;Orig_CCA2_mob=CCA2_no_sub_orig_mob;CCA1_mob=CCA1_no_sub_mob;%Inc
load CCA1_mob_strmANDsbsrfTishrey5787.txt;load CCA2_mob_noStrmNoSbsrfce1.txt;load CCA2_mob_noStrmNoSbsrfce_orig1.txt;x1=CCA1_mob_strmANDsbsrfTishrey5787;x2=CCA2_mob_noStrmNoSbsrfce1;Orig_CCA2_mob=CCA2_mob_noStrmNoSbsrfce_orig1;%Including subsurface, not stream. 
%This is best set
%load CCA1_mob_noStrmWthSbsrfce.txt;load CCA2_mob_noStrmNoSbsrfce1.txt;load CCA2_mob_noStrmNoSbsrfce_orig1.txt;x1=CCA1_mob_noStrmWthSbsrfce;x2=CCA2_mob_noStrmNoSbsrfce1;Orig_CCA2_mob=CCA2_mob_noStrmNoSbsrfce_orig1;
%%%%%%
%The best set Kislev 5786
%load CCA1_mob_noStrmNoSbsrfce1.txt;load CCA2_mob_noStrmNoSbsrfce1.txt;load CCA2_mob_noStrmNoSbsrfce_orig1.txt;x1=CCA1_mob_noStrmNoSbsrfce1;x2=CCA2_mob_noStrmNoSbsrfce1;Orig_CCA2_mob=CCA2_mob_noStrmNoSbsrfce_orig1;%Inc
%%%%Av 5786%%%%%%

%load CCA1_Av_5786_NoStreamNoSubsurface.txt;load CCA2_Elul.txt;load CCA2_El_orig.txt; x1=CCA1_Av_5786_NoStreamNoSubsurface;x2=CCA2_Elul;Orig_CCA2_mob=CCA2_El_orig;%Av_5786_NoStreamWthSubsurface_d=dominant flowpath for current dataset. load CCA2_El_5786_NoStreamNoSubsurface.txt;
%load Xran_mob; x1=Xrand;
%Av_5786_NoStreamWthSubsurface
%load XandY.mat;x1=X;x2=Y;
%%%%%%%
%Test for rand mat
% x1 = load('CCA1_Av_random_rowperm.txt');
% x2 = load('CCA2_Av_5786_NoStreamNoSubsurface.txt');
% Orig_CCA2_mob = load('CCA2_Av_orig.txt');


%PAY ATTENTION: NaNs in CCA1_mob were replaced to 0.0001 instead of NaNs.
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
% %% ========= PCA =========
% %Standardize compounds (feature scaling
% %filename='loc.txt';%Full locations (107) (works also for pca)
% filename='loc_flwpth.txt';%Full locations (107) (works also for pca)
% %filename='GW_name.txt';
% T=readtable(filename);
% varNamesCCA1=T{:,1};
% [x1, coeff1, expl1, keepIdx1] = pca_for_kcca(CCA1_mob, varNamesCCA1, 90, 20);
%%
%A=CCA1_shrt(:,1);
% shuffled_A = A(randperm(length(A)));
% CCA1_shrt(:,1)=shuffled_A;
%filename='loc_name.txt';
% filename='loc_nameNoStream.txt';
%filename='nameShrt.txt';%No manholes
%filename='compounds.txt';%23 compounds
filename='42compounds.txt';
%filename='names_full_mob.txt';%Full compounds-48 (works also for pca)
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
x1 = x1 - mean(x1,'omitnan');
x2 = x2 - mean(x2,'omitnan');

% Normalize variance
 x1 = x1 * sqrt(diag(1./diag(x1'*x1)));
% x2 = x2 * sqrt(diag(1./diag(x2'*x2)));


%%
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
% === Clustering and Visualization ===
Z = [y1, y2];  % 2D embedding

% Plot labeled groups
% figure;
 %gscatter(Z(:,1), Z(:,4), labels);% Z(:,[1 4]) because numeign=3 means three first are by x1 and three last, i.e. 4:6, are by x2. Best correlated are the first ones in each set (i.e.columns 1 and 4).
% % (Z(:,1), Z(:,2), labels);Scatter is based on the vectors in the kernel space that are most correlated. 
% xlabel('KCCA 1'); ylabel('KCCA 2');
% title('KCCA Projection by Sample Type');
% 
%R=corr(x1,y1);%This function checks how the x1 original data are correlated with the latent space of y1.
% %i.e. how they most affect the y1 vector. Indeed chech which pesticides
% %most effect the y1 in the latent space, thus the KCCA
% % Optional: run unsupervised clustering (k-means)
flow3 = categorical(Orig_CCA2_mob(:,1), [1 2 3], {'Runoff','subsurface','Groundwater'});
%flow3 = categorical(Orig_CCA2_mob(:,3), [1 2 3 4 5 6], {'Runoff','subsurface','Groundwater'});
figure(1);
gscatter(Z(:,1), Z(:,4), flow3);  % color by 3-class flowpath
xlabel('KCCA 1'); ylabel('KCCA 2'); title('KCCA Projection by flowpath');

% optional: add compound text labels without affecting colors
hold on
for i = 1:numel(name)
    text(Z(i,1), Z(i,4), ['  ' char(name(i))], 'FontSize',8);
end


%% ============================================================
% === VIOLIN PLOTS: FLOWPATH EFFECT ON COMPOUND KCCA SCORES ===
% =============================================================

%% === Setup for flowpath categories ===
% General (3): 1=Runoff, 2=Subsurface, 3=Groundwater
% codes_general = Orig_CCA2_mob(:,1);
% names_general = {'Runoff','Subsurface','Groundwater'};
% groups_general = categorical(codes_general,[1 2 3],names_general);

% Detailed (6): 1=GW, 2=Pipe, 3=Manhole, 4=Overland, 5=Secondary trench, 6=Primary trench
codes_detailed = Orig_CCA2_mob(:,1);
%names_detailed = {'Groundwater','Subsurface','Overland flow','Channels'};kk=[1 2 4 5];gNames= {"GW", "SB", "OF", "CH"};
%names_detailed = {'Groundwater','Pipe','Manhole','Overland flow','Secondary channels','Primary channels'};kk=[1 2 3 4 5 6];;gNames= {"GW", "P","MN", "OF", "SC", "PC"};
names_detailed = {'Groundwater','Subsurface','Overland flow','Secondary channels','Primary channels'};kk=[1 2 3 4 5];;gNames= {"GW", "SB", "OF", "ST", "PT"};
groups_detailed = categorical(codes_detailed,kk,names_detailed);

% %% === Setup 1: Runoff / Subsurface / Groundwater ===
% figure(2)
% scores1 = y1(:,1);                 
% groups1 = CCA25_mob_orig(:,1);           
% groupNames1 = {'Runoff','subsurface','Groundwater'};
% 
% figure('Units','normalized','Position',[0.3 0.3 0.4 0.5]);
% violinplot(scores1, groups1);
% title('Y₁ distribution per flowpath');
% xlabel('Dominant flowpath');
% ylabel('Canonical variate (Y₁)');
% grid on
% xticks(1:numel(groupNames1));
% xticklabels(groupNames1);
% xtickangle(20);
% 
% % Annotate sample counts
% hold on
% [grpIdx,~] = findgroups(groups1);
% nPer = splitapply(@numel, scores1, grpIdx);
% yl = ylim;
% ytxt = yl(2) - 0.05*range(yl);
% for i = 1:numel(nPer)
%     text(i, ytxt, sprintf('n=%d', nPer(i)), 'HorizontalAlignment','center');
% end
% 
% % Kruskal–Wallis test
% [p_kw1, tbl1, stats1] = kruskalwallis(scores1, groups1, 'off');
% H1 = cell2mat(tbl1(2,5)); % extract H statistic
% text(mean(xlim), yl(2)*0.95, sprintf('KW p = %.2g', p_kw1), ...
%      'HorizontalAlignment','center','FontWeight','bold');
% disp('==== Kruskal–Wallis: Flowpath ====');
% fprintf('Global H = %.2f, p = %.3g\n', H1, p_kw1);
% 
% % Dunn–Šidák post-hoc
% res1 = multcompare(stats1, 'CType','dunn-sidak', 'Display','off');
% 
% % Determine actual internal group order
% u1 = unique(groups1,'stable');
% labels1 = groupNames1(1:numel(u1));  % adapt dynamically
% 
% disp('--- Dunn–Šidák post-hoc results (Flowpath) ---');
% T1 = table(labels1(res1(:,1))', labels1(res1(:,2))', ...
%            res1(:,4), res1(:,6), ...
%            'VariableNames', {'Group1','Group2','Difference','pValue'});
% disp(T1);
% 
% % Heatmap
% k1 = numel(u1);
% P1 = nan(k1);
% for r = 1:size(res1,1)
%     i = res1(r,1); j = res1(r,2);
%     P1(i,j) = res1(r,6); P1(j,i) = res1(r,6);
% end
% figure('Name','Post-hoc map — Flowpath','Units','normalized','Position',[0.2 0.2 0.35 0.5]);
% imagesc(-log10(P1)); colorbar;
% title('-log_{10}(p) pairwise significance');
% xticks(1:k1); yticks(1:k1);
% xticklabels(labels1); yticklabels(labels1);
% xtickangle(30);

%% === Setup 2: Groundwater / Pipe / Manhole / Overland flow / Secondary trench / Primary trench ===
scores2 = y1(:,1);                 
groups2 = codes_detailed;%(:,1)           
%groupNames2 = {'Groundwater','Pipe','Manhole','Overland flow','Secondary trench','Primary trench'};
groupNames2 = names_detailed;%{'Groundwater','Subsurface','Soil','Overland flow','Channels'};

figure('name','many flowpath','Units','normalized','Position',[0.3 0.3 0.4 0.5]);
violinplot(scores2, groups2);
title('Y₁ distribution per sub-flowpath');
xlabel('Sub-flowpath');
ylabel('Canonical variate (Y₁)');
grid on
xticks(1:numel(groupNames2));
xticklabels(groupNames2);
xtickangle(20);

% Annotate sample counts
hold on
[grpIdx,~] = findgroups(groups2);
nPer = splitapply(@numel, scores2, grpIdx);
yl = ylim;
ytxt = yl(2) - 0.05*range(yl);
for i = 1:numel(nPer)
    text(i, ytxt, sprintf('n=%d', nPer(i)), 'HorizontalAlignment','center');
end

% Kruskal–Wallis test
[p_kw_Y1, tbl2, stats_Y1] = kruskalwallis(scores2, groups2, 'off');
%H_Y1 = tbl_Y1{2,5}; 
H_Y1 = cell2mat(tbl2(2,5)); %H2 extract H statistic
text(mean(xlim), yl(2)*0.95, sprintf('KW p = %.2g', p_kw_Y1), ...
     'HorizontalAlignment','center','FontWeight','bold');
disp('==== Kruskal–Wallis: Sub-flowpath ====');
fprintf('Global H = %.2f, p = %.3g\n', H_Y1, p_kw_Y1);

% Dunn–Šidák post-hoc
res2 = multcompare(stats_Y1, 'CType','dunn-sidak', 'Display','off');
u2 = unique(groups2,'stable');
labels2 = groupNames2(1:numel(u2));  % adapt dynamically

disp('--- Dunn–Šidák post-hoc results (Sub-flowpath) ---');
T2 = table(labels2(res2(:,1))', labels2(res2(:,2))', ...
           res2(:,4), res2(:,6), ...
           'VariableNames', {'Group1','Group2','Difference','pValue'});
disp(T2);

% Heatmap
k2 = numel(u2);
P2 = nan(k2);
for r = 1:size(res2,1)
    i = res2(r,1); j = res2(r,2);
    P2(i,j) = res2(r,6); P2(j,i) = res2(r,6);
end
figure('Name','Post-hoc map — Sub-flowpath','Units','normalized','Position',[0.2 0.2 0.35 0.5]);
imagesc(-log10(P2)); 
%colormap(viridis);
colormap(parula);%This enable discriminate colors. If p value is 10^-3 it gives blue, interprating it as low, although it is significant
colorbar;
%caxis([0 5]);%it changes only the colorbar scale due to the power, not the reals values, that significant changes would be perceptional
title('-log_{10}(p) pairwise significance');

xticks(1:k2); yticks(1:k2);
xticklabels(labels2); yticklabels(labels2);
xtickangle(30);
M = -log10(P2);%P2%
% 
% % ============================
% % Add numeric labels
% % ============================
% for i = 1:k2
%     for j = 1:k2
%         if ~isnan(M(i,j))
%             val = M(i,j);
%             if val < 10
%                 txt = sprintf('%.2f', val);   % format: 2 decimals
%             else
%                 txt = sprintf('%.1f', val);   % avoid large clutter
%             end
%             text(j, i, txt, ...
%                 'HorizontalAlignment', 'center', ...
%                 'FontSize', 11, ...
%                 'Color', 'w', ...
%                 'FontWeight','bold');
%         end
%     end
% end
% 
% hold off;
%% ============================================
% === HELLINGER DISTANCE BETWEEN FLOWPATHS ===
% ============================================

% Use the dominant flowpath groups you already defined
groupCodes = codes_detailed; %Orig_CCA2_mob(:,1);     % 1–3 or 1–6 depending on dataset
uniqueGroups = unique(groupCodes,'stable');
%gNames = names_general;              % if using 3 groups
gNames = gNames;           % if using 6 groups – choose manually

% We use the canonical variate Y1(:,1)
vals = y1(:,1);

% Compute group-wise histograms with same binning
nbins = 20;   % adjust if needed; using equal bins ensures comparable densities
edges = linspace(min(vals), max(vals), nbins+1);

% Compute PDFs per group
k = numel(uniqueGroups);
pdfs = zeros(k, nbins);

for i = 1:k
    grp = uniqueGroups(i);
    these = vals(groupCodes == grp);
    pdfs(i,:) = histcounts(these, edges, 'Normalization','probability');
end

% Compute pairwise Hellinger distances
Hmat = zeros(k);

for i = 1:k
    for j = i+1:k
        p = pdfs(i,:);
        q = pdfs(j,:);
        H = sqrt( 1 - sum( sqrt(p.*q) ) );    % standard Hellinger distance
        Hmat(i,j) = H;
        Hmat(j,i) = H;
    end
end




%% === Plot Hellinger heatmap ===
figure('Name','Hellinger distance between flowpath distributions',...
    'Units','normalized','Position',[0.25 0.25 0.35 0.45]);
imagesc(Hmat);
colorbar;
title('Hellinger distance between flowpath Y_1 distributions');
axis square;
xticks(1:k); yticks(1:k);
xticklabels(gNames(1:k)); yticklabels(gNames(1:k));
xtickangle(25);
set(gca,'FontSize',11);


%% ================================================================
% === FLOWPATH CLUSTERS: Do they have distinct Koc by Hellinger? ==
% ================================================================

% --- Use the SAME flowpath grouping you already used ---
groupCodes = codes_detailed; %Orig_CCA2_mob(:,1);     % 1..3 (Runoff/Subsurface/GW) or your actual coding
%gNames     = {"GW", "SB", "OF", "ST", "PT"};% gw=groundwater, sb-subsurface, of-overland flow (field runoff), sc-secondary channel, pc-primary channel

% --- Pull property: Koc (user: "third column") ---
% If in your files Koc is actually column 3, just set idxKoc = 3 or the desired mobility property (Log Kow etc.).
idxKoc = 2; %1-Log Kow, 2-Koc, 3-Solubility 4-Degradability (5-application time).                      
idxKoc=idxKoc+1;
Koc = Orig_CCA2_mob(:, idxKoc);%Koc is general name for mobility properties

% --- Keep only rows where both group and Koc are valid ---
valid = ~isnan(groupCodes) & ~isnan(Koc);
groupCodes = groupCodes(valid);
Koc    = Koc(valid);

% --- Prepare common binning for Koc across groups ---
nbins = 20;  % you can tweak
lo = min(Koc);
hi = max(Koc);
if lo == hi
    warning('All Koc are identical; Hellinger distances will be zero.');
    hi = lo + eps;
end
edges = linspace(lo, hi, nbins+1);

% --- Build PDFs per flowpath cluster (same clusters as your violins) ---
uG = unique(groupCodes,'stable');
k  = numel(uG);
pdfs = zeros(k, nbins);
for i = 1:k
    xi = Koc(groupCodes == uG(i));
    pdfs(i,:) = histcounts(xi, edges, 'Normalization','probability');%Creating the Koc distribution of a group in histogram
end
% ========================================
% (Optional) Violin + Hellinger + KW/Dunn
% ========================================

% --- Pairwise Hellinger distances between flowpath Koc distributions ---
HD_Koc = zeros(k);
for i = 1:k
    for j = i+1:k
        p = pdfs(i,:); %#ok<*NBRAK>  % (MATLAB will accept p = pdfs(i,:);)
        q = pdfs(j,:);
        p = pdfs(i,:);
        q = pdfs(j,:);
        H = sqrt( 1 - sum( sqrt(p .* q) ) );  % 0..1
        HD_Koc(i,j) = H;
        HD_Koc(j,i) = H;
    end
end

% --- Plot Hellinger heatmap: distinctness of Koc by flowpath cluster ---
figure('Name','Hellinger distance of Koc across flowpath clusters',...
       'Units','normalized','Position',[0.26 0.26 0.36 0.45]);
imagesc(HD_Koc);
axis square; colorbar;
title('Hellinger distance: Koc distributions by flowpath');
xticks(1:k); yticks(1:k);
xticklabels(gNames(1:k)); yticklabels(gNames(1:k));
xtickangle(25);
set(gca,'FontSize',11);


doExtras = true;  % set false to skip
if doExtras
    % Violin of Koc by flowpath (same groups)
    figure('Name','Koc by flowpath (violins)', 'Units','normalized','Position',[0.15 0.2 0.45 0.55]);
    violinplot(Koc, groupCodes);
    ylabel('Koc'); xlabel('Flowpath');
    title('Koc distribution per flowpath (same clusters as y_1 violins)');
    grid on; xticks(1:k); xticklabels(gNames(1:k)); xtickangle(20);

    % Annotate n per group
    hold on;
    [grpIdx,~] = findgroups(groupCodes);
    nPer = splitapply(@numel, Koc, grpIdx);
    yl = ylim; ytxt = yl(2) - 0.06*range(yl);
    for i = 1:k, text(i, ytxt, sprintf('n=%d', nPer(i)), 'HorizontalAlignment','center'); end
%=========================================
%     % Kruskal–Wallis + Dunn–Šidák
%=========================================
    [p_kw_Koc, tbl_Koc, stats_Koc] = kruskalwallis(Koc, groupCodes, 'off');
    H_Koc = tbl_Koc{2,5};    
    %H = cell2mat(tbl(2,5));
    text(mean(xlim), yl(2) - 0.12*range(yl), sprintf('KW: H=%.2f, p=%.2g', H, p_kw_Koc), ...
         'FontWeight','bold','HorizontalAlignment','center');

    R = multcompare(stats_Koc, 'CType','dunn-sidak', 'Display','off'); % [i j L diff U p]
    % Make a compact -log10(p) map for the pairwise test
    P = nan(k); 
    for r = 1:size(R,1)
        i = R(r,1); j = R(r,2);
        P(i,j) = R(r,6); P(j,i) = R(r,6);
    end
    %Dunn-Sidak bu mobility
    figure('Name','Post-hoc map: Koc by flowpath','Units','normalized','Position',[0.2 0.2 0.35 0.45]);
    imagesc(-log10(P)); colorbar;
    title('-log_{10}(p) Dunn–Šidák (Koc by flowpath)');
    axis square; xticks(1:k); yticks(1:k);
    xticklabels(gNames(1:k)); yticklabels(gNames(1:k));
    xtickangle(30);
end


%% ============================================================
% === PERMUTATION TEST ADDED ONLY — EXISTING SCRIPT UNCHANGED ===
% ============================================================
% Tests whether the selected mobility property differs among the
% existing dominant-flowpath groups more than expected by chance.
% Existing group sizes are preserved by permuting the group labels.

nPerm_mob = 5000;
rng(1,'twister');

% Observed KW H for the selected mobility property
[p_kw_perm_reference, tbl_kw_perm_reference] = kruskalwallis(Koc, groupCodes, 'off');
H_Koc_perm_reference = tbl_kw_perm_reference{2,5};

% Observed global mean pairwise Hellinger distance
upperMask_HD = triu(true(k),1);
meanHD_obs = mean(HD_Koc(upperMask_HD),'omitnan');

% Null distributions
H_kw_null = nan(nPerm_mob,1);
meanHD_null = nan(nPerm_mob,1);
maxHD_null = nan(nPerm_mob,1);
HD_null = nan(k,k,nPerm_mob);

for ip = 1:nPerm_mob

    % Permute existing flowpath labels; group sizes therefore stay identical
    permGroups = groupCodes(randperm(numel(groupCodes)));

    % ---- Permutation KW H ----
    [~, tbl_perm] = kruskalwallis(Koc, permGroups, 'off');
    H_kw_null(ip) = tbl_perm{2,5};

    % ---- Permutation Hellinger distances ----
    pdfs_perm = zeros(k,nbins);
    for ig = 1:k
        xi_perm = Koc(permGroups == uG(ig));
        pdfs_perm(ig,:) = histcounts(xi_perm, edges, ...
            'Normalization','probability');
    end

    HDp = zeros(k);
    for ii = 1:k
        for jj = ii+1:k
            pp = pdfs_perm(ii,:);
            qq = pdfs_perm(jj,:);
            hp = sqrt(1 - sum(sqrt(pp .* qq)));
            HDp(ii,jj) = hp;
            HDp(jj,ii) = hp;
        end
    end

    HD_null(:,:,ip) = HDp;

    tmpHD = HDp(upperMask_HD);
    meanHD_null(ip) = mean(tmpHD,'omitnan');
    maxHD_null(ip) = max(tmpHD);
end

% Global permutation p-values (+1 correction)
p_kw_permutation = (1 + sum(H_kw_null >= H_Koc_perm_reference)) ...
                   / (nPerm_mob + 1);

p_meanHD_permutation = (1 + sum(meanHD_null >= meanHD_obs)) ...
                       / (nPerm_mob + 1);

% Pairwise HD permutation p-values and max-statistic FWER
P_HD_perm = nan(k);
P_HD_FWER = nan(k);

for ii = 1:k
    for jj = ii+1:k

        obsHD = HD_Koc(ii,jj);
        pairNull = squeeze(HD_null(ii,jj,:));

        pRaw = (1 + sum(pairNull >= obsHD)) / (nPerm_mob + 1);
        pFWER = (1 + sum(maxHD_null >= obsHD)) / (nPerm_mob + 1);

        P_HD_perm(ii,jj) = pRaw;
        P_HD_perm(jj,ii) = pRaw;

        P_HD_FWER(ii,jj) = pFWER;
        P_HD_FWER(jj,ii) = pFWER;
    end
end

fprintf('\n==== PERMUTATION TEST: mobility by flowpath ====\n');
fprintf('Permutations = %d\n',nPerm_mob);
fprintf('KW observed H = %.4f, permutation p = %.6g\n', ...
    H_Koc_perm_reference,p_kw_permutation);
fprintf('Mean HD observed = %.4f, global permutation p = %.6g\n', ...
    meanHD_obs,p_meanHD_permutation);

% Pairwise results table
nHDpairs = k*(k-1)/2;
HDPermutationTable = table('Size',[nHDpairs 5], ...
    'VariableTypes',{'string','string','double','double','double'}, ...
    'VariableNames',{'Group1','Group2','HD','PermutationP','FWER_P'});

rr = 0;
for ii = 1:k
    for jj = ii+1:k
        rr = rr + 1;
        HDPermutationTable.Group1(rr) = string(gNames{ii});
        HDPermutationTable.Group2(rr) = string(gNames{jj});
        HDPermutationTable.HD(rr) = HD_Koc(ii,jj);
        HDPermutationTable.PermutationP(rr) = P_HD_perm(ii,jj);
        HDPermutationTable.FWER_P(rr) = P_HD_FWER(ii,jj);
    end
end

disp(HDPermutationTable);

% Raw pairwise HD permutation p-value heatmap
Praw_plot = nan(k);
mask = isfinite(P_HD_perm);
Praw_plot(mask) = -log10(max(P_HD_perm(mask),realmin));

figure('Name','HD permutation p: mobility by flowpath', ...
       'Units','normalized','Position',[0.20 0.20 0.38 0.46]);
h = imagesc(Praw_plot);
set(h,'AlphaData',isfinite(Praw_plot));
set(gca,'Color','k');
axis square;
colormap(parula);
colorbar;
title('-log_{10}(p) HD permutation: mobility by flowpath');
xticks(1:k); yticks(1:k);
xticklabels(gNames(1:k)); yticklabels(gNames(1:k));
xtickangle(25);

% Max-statistic FWER heatmap
Pfwer_plot = nan(k);
mask = isfinite(P_HD_FWER);
Pfwer_plot(mask) = -log10(max(P_HD_FWER(mask),realmin));

figure('Name','HD max-stat FWER: mobility by flowpath', ...
       'Units','normalized','Position',[0.22 0.22 0.38 0.46]);
h = imagesc(Pfwer_plot);
set(h,'AlphaData',isfinite(Pfwer_plot));
set(gca,'Color','k');
axis square;
colormap(parula);
colorbar;
title('-log_{10}(FWER p) HD permutation: mobility by flowpath');
xticks(1:k); yticks(1:k);
xticklabels(gNames(1:k)); yticklabels(gNames(1:k));
xtickangle(25);


nPairs = size(R,1);
DunnTable = table('Size',[nPairs 6], ...
    'VariableTypes', {'string','string','double','double','double','double'}, ...
    'VariableNames', {'Group_i','Group_j','DiffMeanRank','LowerCI','UpperCI','pValue'});

% --- Dunn–Šidák post-hoc results (Flowpath groups) ---
R = multcompare(stats_Koc, 'CType','dunn-sidak', 'Display','off');

u = unique(groupCodes,'stable');
labels = names_detailed(1:numel(u));   % dynamic group-name mapping

% disp('--- Dunn–Šidák post-hoc results (Flowpath) ---');
% T = table(labels(R(:,1))', labels(R(:,2))', ...
%     R(:,4), R(:,6), ...
%     'VariableNames', {'Group1','Group2','Difference','pValue'});
% disp(T);
% %=============================
% % Add numeric labels
% % ============================

% Label each cell with numeric values

% Label each cell with numeric values


% hold on;
% P1=-log10(P)%P%
% % === Add numeric labels to each cell ===
% for ii = 1:k
%     for jj = 1:k
%         if ~isnan(P1(ii,jj))
%             text(jj, ii, sprintf('%.2f', P1(ii,jj)), ...
%                 'HorizontalAlignment','center', ...
%                 'FontWeight','bold', ...
%                 'Color','w', ...
%                 'FontSize',11);
%         end
%     end
% end
% 
% hold off;
%%
FindStatistic
    %% regression for other variables

    controlNames = {'Log Kow','Koc','Solubility','Degradability'};%};
figure('Units','normalized','Position',[0.05 0.15 0.9 0.75]);
for i = 1:length(controlNames)
    subplot(2,2,i);
    scatter(x2(:,i+1), y1(:,1), 40, 'filled', 'MarkerFaceAlpha',0.6);
    hold on;
[rho(i),pval(i)] = corr(x2(:,i+1), y1(:,1), 'type','Spearman');
    % Add linear fit
    h = lsline;
    h.Color = [0 0 0]; % black line
    h.LineWidth = 1.5;
R2(i) = corr(x2(:,i+3), y1(:,1), 'type','Pearson');
    % Axis labels and title
    xlabel(controlNames{i}, 'Interpreter','none');
    ylabel('Canonical Component 1');
    %title(['\rho = ' num2str(rho(i), '%.2f') ',  p = ' num2str(pval(i), '%.3f')]);
    title(['R^2 = ' num2str(R2(i), '%.2f') ',  p = ' num2str(pval(i), '%.3f')]);


    % Optional: beautify
    box on; grid on;
    set(gca,'FontSize',11);
    hold off;
end

%% Correlation table
controlNames = {'Koc','Log Kow','Solubility','Degradability'};%,'Application_time'};

% Create logical column marking significance (p < 0.05)
Significant = pval < 0.05;

% Create table
CorrTable = table(controlNames', rho',R2', pval', Significant', ...
    'VariableNames', {'Variable','Rho', 'R^2','P_value','Significant'});



%% === Inputs assumed already defined above ===
% y1  : [n x k] canonical variates (use y1(:,1))
% x2  : [n x 3] integer labels: [Flowpath(1..7), Timing(1..4), Field(1..2)]
% scores   = y1(:,1);
% Koc          = x2(:,5);
% log_Kow       = x2(:,4);
% Solubility   = x2(:,6);
% Degradability= x2(:,7);
% ApplicationTime= x2(:,8);


%% ===== Helpers =====
function add_counts(ax, vals, groups)
    axes(ax); hold on
    [grpIdx,~] = findgroups(groups);
    nPer = splitapply(@numel, vals, grpIdx);
    yl = ylim(ax); ytxt = yl(2) - 0.06*range(yl);
    for i = 1:numel(nPer)
        text(ax, i, ytxt, sprintf('n=%d', nPer(i)), 'HorizontalAlignment','center');
    end
    hold off
end

function annotate_kw(ax, vals, ids, labelStr)
    [p,~,~] = kruskalwallis(vals, ids, 'off');              % global KW p
    xl = xlim(ax);  yl = ylim(ax);
    xpos = mean(xl);
    ypos = yl(2) - 0.10*range(yl);
    text(ax, xpos, ypos, sprintf('%s: KW p = %.2g', labelStr, p), ...
         'FontWeight','bold','HorizontalAlignment','center');

 [p_kw,~,stats] = kruskalwallis(vals, ids, 'off');
    fprintf('Global KW p = %.3g\n', p_kw);
    R = multcompare(stats, 'CType','dunn-sidak', 'Display','off'); % [i j L diff U p]
    T = table( names(R(:,1))', names(R(:,2))', R(:,4), R(:,6), ...
               'VariableNames', {'Group1','Group2','Difference','pValue'});
    disp(T);
end

% function kw_posthoc(vals, ids, names)
%     [p_kw,~,stats] = kruskalwallis(vals, ids, 'off');
%     fprintf('Global KW p = %.3g\n', p_kw);
%     R = multcompare(stats, 'CType','dunn-sidak', 'Display','off'); % [i j L diff U p]
%     T = table( names(R(:,1))', names(R(:,2))', R(:,4), R(:,6), ...
%                'VariableNames', {'Group1','Group2','Difference','pValue'});
%     disp(T);
% end
