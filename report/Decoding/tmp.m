clear; close all;clc;

%% load Data
root = '/Users/eyalsoreq/Dropbox/PhD/CURRENT/papers/WM/wm-decoding-project/results';
outDir = [root,'/report/Decoding'];
mkdir(outDir);
load([root,'/Data/ModelPerf', '/SVMECOC.mat']);
setenv('PATH', [getenv('PATH') ':/usr/local/bin'] )
format shortG
labels = {'Null','oob','cv','test'};
factors = unique(T_factors.factor);
metric =  unique(T_factors.metric);
comp = unique(T_factors.set)';
comp = comp([3,1,2]);
compL = {'WB','DG','SG'};


%% plot boxplots 
cmap = hsv(12);
for m= metric'
    ix = c3nl.strDetect(T_factors.factor,'state')&c3nl.strDetect(T_factors.metric,m{1});
    iy = c3nl.strDetect(T_factors.factor,'domain')&c3nl.strDetect(T_factors.metric,m{1});
    xx = T_factors(ix,:);
    yy = T_factors(iy,:);
    [y,x,h,w] = get.fancyGrid(ones(numel(comp),1),0.02,0.2,'stretch');
    figure(1);clf;
    ax = [];
    
    for jj= 1:numel(comp)
        ax.(['p' num2str(jj)]) = axes('position',[x(jj),y(jj),w,h]);
        hold ( ax.(['p' num2str(jj)]),'on');
        ix = c3nl.strDetect(xx.set,comp{jj});
        for ii=1:4
            plot(xx{ix,ii+4},yy{ix,ii+4},'Marker','o','LineStyle','none','MarkerFaceColor',min(1,cmap(ii,:)+0.1).^.3,'MarkerEdgeColor','none','MarkerSize',4);       
        end
        for ii=1:4
            ix = c3nl.strDetect(xx.set,comp{jj});
            plot(nanmean(xx{ix,ii+4}),nanmean(yy{ix,ii+4}),'Marker','o','LineStyle','none','MarkerFaceColor',cmap(ii,:),'MarkerEdgeColor','k','MarkerSize',12,'Parent',ax.(['p' num2str(jj)]));       
        end
        for ii=10:10:90
           line([0,100],[ii,ii],'linewidth',0.5,'color',[0.7,0.7,0.7]);
           line([ii,ii],[0,100],'linewidth',0.5,'color',[0.7,0.7,0.7]);
        end
        grid off
        ax.(['p' num2str(jj)]).Color= [0.9,0.9,0.9]; 
        xlabel('Stage F1Score')
        if strcmpi(compL{jj},'WB');ylabel('Domain F1Score');
        else
            ax.(['p' num2str(jj)]).YTickLabel =''; 
        end
        axis ([0,100,0,100])
        title(compL{jj})
    end
    save.pdf(sprintf('%s%sscatter_%s.ai',outDir ,filesep,m{1}),30,11)
end

%% extract empirical p-values for classification vs null 

T_factors.id = repmat(1:1000,1,18)';
T_factors.set = categorical(T_factors.set,unique(T_factors.set),{'DG','SG','WB'});
T_factors.factor = categorical(T_factors.factor,unique(T_factors.factor),{'Domain','Load','Stage'});%
factors = unique(T_factors.factor);
T_factors.metric = categorical(T_factors.metric);
Tg = table();
for f = factors'
    for m = metric'
        for s = unique(T_factors.set)'
            ix = T_factors.metric==m{1}&T_factors.factor==f&T_factors.set==s;
            [p,SEM,ts,CI] =get.empiricalP(T_factors.Null(ix),T_factors.replication(ix),19);
            [D,sd,m1,m2,sd1,sd2] = get.mannwhitneyR(T_factors.replication(ix),T_factors.Null(ix));
            Tg = [Tg;table(m,f,s,m1,sd1,D,p,CI)];
        end
    end
end
Tg.Properties.VariableNames = {'Metric','Factor','ROI_set','mean_F1', 'std_F1','MannWhitneyR', 'p_value','CI'};
caption = ['classification effect size compared to the permuted null effect over metrics (activation and connectivity)',...
           ' Factor and set. Empirical p is calculated as the number of null models that are better than the mean/min F1 of the true model ($p_mean,p_min$).']; 

plot.mdl2Table(Tg(Tg.Factor~='Load',:),[outDir,filesep,'ClassificationEffect.tex'],'latex',caption,'ClassificationEffect')

Treport = grpstats(T_factors,{'set','factor','metric'},{'mean'},'DataVars',{'Null','oob','test','replication'});
Treport.Properties.VariableNames = {'set','factor','metric','GroupCount','Null','oob','CV','Test'};
plot.mdl2Table(Treport(Treport.factor~='Load',:),[outDir,filesep,'ClassificationMeanPerf.tex'],'latex','Mean performance of different classification models','ClassificationMeanPerf')

%% Comparing factors and ROI using non-parametric kruskal wallis testing 

Tg = grpstats(T_factors,{'set','factor','metric','perm'},{'mean'},'DataVars',{'Null','replication'});
factors = categories(T_factors.factor);
factors(2)=[];
for f = factors'
    tmp = Tg(Tg.factor==f{1},:);% extract the factor perfomance 
    tmp.diff = tmp.mean_replication-tmp.mean_Null; % calculate the pairwise distnace from null
    [p,tbl,stats] = kruskalwallis(tmp.diff,tmp.metric,'on');
    caption = sprintf(['Kruskal Wallis non-parametric (distribution free) test comparing '...
        'classification accuracy distance from null between activation and connectivity across sets for decoding %s '...
        '(meanRank$_{%s}$=%4.2f, meanRank$_{%s}$=%4.2f,meanRank$_{%s}$=%4.2f)'],...
        f{1},'Activity',stats.meanranks(1),'Connectivity',stats.meanranks(2));
    tbl = c3nl.roundTable(cell2table(tbl(2:end,:),'VariableNames',matlab.lang.makeValidName(tbl(1,:))),4);
    tbl.etaSquared = zeros(height(tbl),1);
    tbl.etaSquared(1) = tbl.SS(1)/tbl.SS(end);
    plot.mdl2Table(tbl,sprintf('%s%s%s_metric_kruskalwallis.tex',outDir,filesep,f{1}),'latex',caption,'null')
    tbl = table();
    for s= unique(tmp.set)'
        tbl1 = get.ranksum(tmp.diff(tmp.set==s),tmp.metric(tmp.set==s));
        tbl1.Set = {char(s)};
        tbl = [tbl;tbl1];
    end
    caption = 'Post-hoc rank-sum multicomparison';
    plot.mdl2Table(c3nl.roundTable(tbl,4),sprintf('%s%s%s_PHROIkruskalwallis.tex',outDir,filesep,f{1}),'latex',caption,'null');
end
%% compare domain models across loads export empirical p-values for distnace from null 
load([root,'/Data/ModelPerf', '/loadPref.mat']);

T_factors.id = repmat(1:1000,1,6)';
T_factors.set = categorical(T_factors.set,unique(T_factors.set),{'DG','SG','WB'});
T_factors.metric = categorical(T_factors.metric);
Tg = table();

for m = unique(T_factors.metric)'
    for s = unique(T_factors.set)'
        ix = T_factors.metric==m&T_factors.set==s;
        [p,SEM,ts,CI] =get.empiricalP(T_factors.Null(ix),T_factors.replication(ix),19);
        [D,sd,m1,m2,sd1,sd2] = get.mannwhitneyR(T_factors.replication(ix),T_factors.Null(ix));
        Tg = [Tg;table(m,f,s,m1,sd1,D,p,CI)];
    end
end

Tg.Properties.VariableNames = {'Metric','Factor','ROI_set','mean_F1', 'std_F1','MannWhitneyR', 'p_value','CI'};
caption = ['classification effect size compared to the permuted null effect over metrics (activation and connectivity)',...
           ' and set. Empirical p is calculated as the number of null models that are better than the mean F1score of the true model.']; 
plot.mdl2Table(Tg,[outDir,filesep,'ClassificationEffectdomainbyload.tex'],'latex',caption,'ClassificationEffectdomainbyload')


Treport = grpstats(T_factors,{'set','metric'},{'mean'},'DataVars',{'Null','oob','test','replication'});
Treport.Properties.VariableNames = {'set','metric','GroupCount','Null','oob','CV','Test'};
plot.mdl2Table(Treport,[outDir,filesep,'ClassificationMeanPerfdomainbyload.tex'],'latex','Mean performance of different classification models','ClassificationMeanPerf')

% Now establish significant distance between loads using non-parametric kruskal wallis testing 
M_perf.set = categorical(M_perf.set,unique(M_perf.set),{'DG','SG','WB'});

Tg = grpstats(M_perf,{'perm','set','metric'},'mean','DataVars',M_perf.Properties.VariableNames(5:end));
measures = categorical(Tg.metric,{'PPI','BA'},{'Connectivity','Activity'});
for f = unique(measures)'
    ix = measures==f;
    X = Tg{ix,6:end}-Tg{ix,5}; % calculate the pairwise distnace from null
    [p,tbl,stats] = kruskalwallis(reshape(X(:,1:3),[],1),categorical(repelem({'low','med','high'},nnz(ix)),{'low','med','high'}),'off');
    caption = sprintf(['Kruskal Wallis non-parametric (distribution free) test comparing '...
    'classification accuracy distance from null across ROI sets for decoding %s across loads '...
    '(meanRank$_{%s}$=%4.2f, meanRank$_{%s}$=%4.2f,meanRank$_{%s}$=%4.2f)'],...
    char(f),stats.gnames{1},stats.meanranks(1),stats.gnames{2},stats.meanranks(2),stats.gnames{3},stats.meanranks(3));
    tbl = c3nl.roundTable(cell2table(tbl(2:end,:),'VariableNames',matlab.lang.makeValidName(tbl(1,:))),4);
    tbl.etaSquared = zeros(height(tbl),1);
    tbl.etaSquared(1) = tbl.SS(1)/tbl.SS(end);    
    plot.mdl2Table(tbl,sprintf('%s%s%s_Study1kruskalwallis.tex',outDir,filesep,char(f)),'latex',caption,'null')
    [c,m,h,gnames] =multcompare(stats,'CType','bonferroni');
    PHtbl = [array2table(categorical(c(:,1:2),1:numel(gnames),gnames),'VariableNames',{'C1','C2'}),array2table(c(:,3:end))];
    PHtbl.Properties.VariableNames = {'Group1','Group2','lower_limit','mean_diff','upper_limit','p_value'};
    caption = 'Post-hoc test bonferroni corrected for multicomparison';
    plot.mdl2Table(c3nl.roundTable(PHtbl,4),sprintf('%s%s%s_PHLoadkruskalwallis.tex',outDir,filesep,char(f)),'latex',caption,'null');
    [p,h,stats] = ranksum(X(:,4),X(:,5));
    [ranks, tieadj] = tiedrank([X(:,4);X(:,5)]);
    ranks = ranks.*[~repelem([0,1],size(X,1),1);repelem([0,1],size(X,1),1)];ranks(ranks==0)=nan;
    meanranks = nanmean(ranks);
    gp = {'Low','High'};
    N = size(X,1);
    tbl = table(gp(1),gp(2),N,meanranks(1),meanranks(2),stats.zval,p,stats.zval/sqrt(N),'VariableNames',{'Group1','Group2','N','rankMean1','rankMean2','Z','p_value','Effect_Size'});
    caption = sprintf('Wilcoxon rank sum non-parametric (distribution free) test comparing low and high classification accuracy distance from null for decoding visual domains using %s (meanRank$_{low}$=%4.2f, meanRank$_{high}$=%4.2f)',char(f),meanranks(1),meanranks(2));
    plot.mdl2Table(c3nl.roundTable(tbl,4),sprintf('%s%s%s_Study2LoadRankSum.tex',outDir,filesep,char(f)),'latex',caption,'null')   
end

sets= unique(M_perf.set);
clear XBA
for ii=1:numel(sets)
ixBA = c3nl.strDetect(M_perf.metric,'BA')&M_perf.set==sets(ii);
ixPPI = c3nl.strDetect(M_perf.metric,'PPI')&M_perf.set==sets(ii);
XBA{ii,:} = {M_perf.Low_cv(ixBA),M_perf.Med_cv(ixBA),M_perf.High_cv(ixBA)};
%XPPI{ii,1} = {M_perf.Low_cv(ixPPI),M_perf.Med_cv(ixPPI),M_perf.High_cv(ixPPI)};
end
cmap = hsv(12);

cmap =[[256,0,0];0,176,254;0,104,253;0,40,168]./256;

plot.histfit({M_perf.Null(ixBA), M_perf.Low_cv(ixBA),M_perf.Med_cv(ixBA),M_perf.High_cv(ixBA)},{{'CV'},{'Low','Med','High'}},[0,100.5],cmap,1,1)
%legend({'null','Low','Med','High'},'Location','northoutside','Orientation','horizontal')
ax = gca;
ax.Color = [0.9,0.9,0.9];
ax.XTickLabel = '';
save.pdf(sprintf('%s%s_histFit_CV_%s.ai',outDir ,filesep,'BA'),14,3)

plot.histfit({M_perf.Null(ixBA), M_perf.Low_R(ixBA),M_perf.High_R(ixBA)},{{'Test'},{'Low','High'}},[0,100.5],cmap([1:2,4],:),1,2)
ax = gca;
ax.Color = [0.9,0.9,0.9];
save.pdf(sprintf('%s%s_histFit_Test_%s.ai',outDir ,filesep,'BA'),14,3)

plot.histfit({M_perf.Null(ixPPI), M_perf.Low_cv(ixPPI),M_perf.Med_cv(ixPPI),M_perf.High_cv(ixPPI)},{{'CV'},{'Low','Med','High'}},[0,100.5],cmap,1,2);
%legend({'null','Low','Med','High'},'Location','northoutside','Orientation','horizontal')
ax = gca;
ax.Color = [0.9,0.9,0.9];
ax.XTickLabel = '';
save.pdf(sprintf('%s%s_histFit_CV_%s.ai',outDir ,filesep,'PPI'),14,3)

plot.histfit({M_perf.Null(ixPPI), M_perf.Low_R(ixPPI),M_perf.High_R(ixPPI)},{{'Test'},{'Low','High'}},[0,100.5],cmap([1:2,4],:),1,2)
ax = gca;
ax.Color = [0.9,0.9,0.9];
save.pdf(sprintf('%s%s_histFit_Test_%s.ai',outDir ,filesep,'PPI'),14,3)


%% Compare domain models across stages produce empirical p-values

T_perf.id = repmat(1:1000,1,6)';
T_perf.set = categorical(T_perf.set,unique(T_perf.set),{'DG','SG','WB'});
T_perf.metric = categorical(T_perf.metric);
Tg = table();

for m = unique(T_perf.metric)'
    for s = unique(T_perf.set)'
        ix = T_perf.metric==m&T_perf.set==s;
        [p,SEM,ts,CI] =get.empiricalP(T_perf.Null(ix),T_perf.replication(ix),19);
        [D,sd,m1,m2,sd1,sd2] = get.mannwhitneyR(T_perf.replication(ix),T_perf.Null(ix));
        Tg = [Tg;table(m,f,s,m1,sd1,D,p,CI)];
    end
end

Tg.Properties.VariableNames = {'Metric','Factor','ROI_set','mean_F1', 'std_F1','MannWhitneyR', 'p_value','CI'};
caption = ['classification effect size compared to the permuted null effect over metrics (activation and connectivity)',...
           ' and set. Empirical p is calculated as the number of null models that are better than the mean/min F1 of the true model ($p_mean,p_min$).']; 
plot.mdl2Table(Tg,[outDir,filesep,'ClassificationEffectdomainbystage.tex'],'latex',caption,'ClassificationEffect')


Treport = grpstats(T_perf,{'set','metric'},{'mean'},'DataVars',{'Null','oob','test','replication'});
Treport.Properties.VariableNames = {'set','metric','GroupCount','Null','oob','CV','Test'};
plot.mdl2Table(Treport,[outDir,filesep,'ClassificationMeanPerfdomainbystage.tex'],'latex','Mean performance of different classification models','ClassificationMeanPerf')




%% domain by stages by sets models 

sets = unique(T_domainByStages.set);
metric = unique(T_domainByStages.metric);
models= unique(T_domainByStages.model);
subgroup= unique(T_domainByStages.subgroup);
DS = zeros(1000,5,3,3,2);
names = cell(5,3,3,2);
for m=1:numel(metric)
    for ii=1:numel(sets)
        for jj=1:numel(models)
            for k=1:numel(subgroup)
                ix = c3nl.strDetect(T_domainByStages.metric,metric{m})&c3nl.strDetect(T_domainByStages.subgroup,subgroup{k})&c3nl.strDetect(T_domainByStages.model,models{jj})&c3nl.strDetect(T_domainByStages.set,sets{ii});
                DS(:,jj,k,ii,m) =  T_domainByStages.F1_R(ix);
            end
        end
    end
end

% plot the stripboards
for m=1:numel(metric)
    for ii=1:numel(sets)
        plot.stripboard(DS(:,:,:,ii,m),{{'E','M','P','G','N'},{'E','M','P'}},1,[0,100],[1,1,0]);
        save.pdf(sprintf('%s/SB_%s_%s.ai',outDir,metric{m},sets{ii}),5,5);
    end
end

%% match missmatch analysis 


tmp = T_domainByStages(~c3nl.strDetect(T_domainByStages.model,'null'),:);
tmp.gp = categorical(tmp.model).*categorical(tmp.subgroup);
tmp = unstack(table(tmp.set,tmp.metric,tmp.perm,tmp.mc,tmp.gp,tmp.F1_R,'VariableNames',{'set','metric','perm','mc','gp','F1'}),'F1',{'gp'},'AggregationFun',@mean);
tmp{:,5:7}= tmp{:,5:7}-tmp{:,14:16};
tmp{:,8:10}= tmp{:,8:10}-tmp{:,14:16};
tmp{:,11:13}= tmp{:,11:13}-tmp{:,14:16};
Match  = mean(tmp{:,[5,9,13]},2);
MissMatch = mean(tmp{:,[6:8,10:12]},2);

sets = unique(tmp.set);
metric = unique(tmp.metric);
cmap = [45,234,121;57,181,74;0,104,56]./256;
c = 1
for m= metric'
    ix = c3nl.strDetect(tmp.metric,m{1});
    figure(c);clf;
    ax = gca;hold all;
    S = categorical(tmp.set);
    for ii=1:numel(sets)
        ix1 = ix &( S == sets(ii));
        plot(Match(ix1),MissMatch(ix1),'Marker','o','LineStyle','none','MarkerFaceColor',min(1,cmap(ii,:)+0.1).^.7,'MarkerEdgeColor','none','MarkerSize',4);       
        hold (ax,'on');
    end
    for ii=1:numel(sets)
        ix1 = ix &( S == sets(ii));
        plot(mean(Match(ix1)),mean(MissMatch(ix1)),'Marker','o','LineStyle','none','MarkerFaceColor',cmap(ii,:),'MarkerEdgeColor','k','MarkerSize',12);       
    end
    c = c+1;
    
    grid off
    ax.Color = [0.9,0.9,0.9];
    for ii=-50:10:20
       line([-50,20],[ii,ii],'linewidth',0.5,'color',[0.7,0.7,0.7]);
       line([ii,ii],[-50,20],'linewidth',0.5,'color',[0.7,0.7,0.7]);
    end
    line([-30,20],[0,0],'linewidth',0.75,'color',zeros(1,3));
    line([0,0],[-50,10],'linewidth',0.75,'color',zeros(1,3));
    axis([-30,20,-50,10])
    %legend(sets)
    save.pdf(sprintf('%s%smissmatch_scatter_%s.ai',outDir ,filesep,m{1}),12,7)
end

%% test significance of metric (connectivity vs activation) 

tmp = T_domainByStages(~c3nl.strDetect(T_domainByStages.model,'null'),:);
tmp.gp = categorical(tmp.model).*categorical(tmp.subgroup);
tmp = unstack(table(tmp.set,tmp.metric,tmp.perm,tmp.mc,tmp.gp,tmp.F1_R,'VariableNames',{'set','metric','perm','mc','gp','F1'}),'F1',{'gp'},'AggregationFun',@mean);
tmp = grpstats(tmp,{'set','metric','perm'},{'mean'},'DataVars',tmp.Properties.VariableNames(5:end));
tmp.set = categorical(tmp.set,unique(tmp.set),{'DG','SG','WB'});

tmp{:,5:7}= tmp{:,5:7}-tmp{:,14:16};
tmp{:,8:10}= tmp{:,8:10}-tmp{:,14:16};
tmp{:,11:13}= tmp{:,11:13}-tmp{:,14:16};
Match  = mean(tmp{:,[5,9,13]},2);
MissMatch = mean(tmp{:,[6:8,10:12]},2);
ix = c3nl.strDetect(tmp.metric,'BA');
[p,tbl,stats] =kruskalwallis(Match,tmp.metric,'off');
caption = sprintf('Comparing activation and connectivity Match effects (meanrank_{activation}=%2.3f,meanrank_{connectivity}=%2.3f)',stats.meanranks(1),stats.meanranks(2))   ; 
tbl = c3nl.roundTable(cell2table(tbl(2:end,:),'VariableNames',matlab.lang.makeValidName(tbl(1,:))),4);
plot.mdl2Table(tbl,[outDir,filesep,'matchMetric.tex'],'latex',caption,'matchMetric');

[p,tbl,stats] =kruskalwallis(MissMatch,tmp.metric,'off');
caption = sprintf('Comparing activation and connectivity MisMatch effects (meanrank_{activation}=%2.3f,meanrank_{connectivity}=%2.3f)',stats.meanranks(1),stats.meanranks(2))   ; 
tbl = c3nl.roundTable(cell2table(tbl(2:end,:),'VariableNames',matlab.lang.makeValidName(tbl(1,:))),4);
plot.mdl2Table(tbl,[outDir,filesep,'mismatchMetric.tex'],'latex',caption,'mismatchMetric');
tbl = [table(repelem({'Connectivity','Activity'},[3,3])','VariableNames',{'Metric'}), [get.signtest(Match(~ix),categorical(tmp.set(~ix))); get.signtest(Match(ix),categorical(tmp.set(ix)))]];
caption = 'Match sign test for zero median of a single population to investigate directionality within each set'  ; 
plot.mdl2Table(tbl,[outDir,filesep,'matchSets.tex'],'latex',caption,'matchSets');
tbl = [table(repelem({'Connectivity','Activity'},[3,3])','VariableNames',{'Metric'}), [get.signtest(MissMatch(~ix),categorical(tmp.set(~ix))); get.signtest(MissMatch(ix),categorical(tmp.set(ix)))]];
caption = 'Mismatch sign test for zero median of a single population to investigate directionality within each set'  ; 
plot.mdl2Table(tbl,[outDir,filesep,'mismatchSets.tex'],'latex',caption,'mismatchSets');



%% test manipulation 

load([root,'/Data/ModelPerf', '/manipPerf.mat']);
plot.bar(M_perf.test,categorical(M_perf.set).*categorical(M_perf.metric),{'manipulation'},jet(6))



labels = {'Null','oob','cv','test'};
factors = unique(M_perf.factor);
metric =  unique(M_perf.metric);
comp = unique(M_perf.set)';
comp = comp([3,1,2]);
compL = {'WB','DG','SG'};
for m= metric'
for f =factors'
ix = c3nl.strDetect(M_perf.factor,f{1})&c3nl.strDetect(M_perf.metric,m{1});
tmp = M_perf(ix,:);
xx = zeros(1000,3,3);% measures , distributions, sets
for ii = 1:numel(comp)
    idx =  c3nl.strDetect(tmp.set,comp{ii});
    xx(:,:,ii) = tmp{idx,5:7};    
end
plot.boxplot(xx,{labels,compL},[],hsv(12),2,3,1,[0,100.5],1,'F1');
save.pdf(sprintf('%s%sbox_%s_%s.ai',outDir ,filesep,m{1},f{1}),8,6)
end
end


M_perf.id = repmat(1:1000,1,6)';
M_perf.set = categorical(M_perf.set,unique(M_perf.set),{'DG','SG','WB'});
M_perf.metric = categorical(M_perf.metric);
set = unique(M_perf.set);
Tg = table();
for f = factors'
    for m = metric'
        for s = set'
            ix = M_perf.metric==m{1}&M_perf.set==s;           
            [p,SEM,ts,CI] =get.empiricalP(M_perf.Null(ix),M_perf.replication(ix),19);
            [D,sd,m1,m2,sd1,sd2] = get.mannwhitneyR(M_perf.replication(ix),M_perf.Null(ix));
            Tg = [Tg;table(m,f,s,m1,sd1,D,p,CI)];
        end
    end
end
Tg.Properties.VariableNames = {'Metric','Factor','ROI_set','mean_F1', 'std_F1','MannWhitneyR', 'p_value','CI'};
caption = ['classification effect size compared to the permuted null effect over metrics (activation and connectivity)',...
           ' Factor and set. Empirical p is calculated as the number of null models that are better than the mean F1 score of the true model.']; 

plot.mdl2Table(Tg,[outDir,filesep,'ClassificationManipEffect.tex'],'latex',caption,'ClassificationManipEffect')


Treport = grpstats(M_perf,{'set','factor','metric'},{'mean'},'DataVars',{'Null','oob','test'});
Treport.Properties.VariableNames = {'set','factor','metric','GroupCount','Null','oob','CV'};
plot.mdl2Table(Treport,[outDir,filesep,'ClassificationManipulationMeanPerf.tex'],'latex','Mean performance of different classification models for the manipulation factor on study 2','ClassificationmaniMeanPerf')


c3nl.copy(c3nl.select('name','*.tex'),'/Users/eyalsoreq/Dropbox/PhD/CURRENT/papers/WM/wm-decoding-project/supplamentryMaterial')
