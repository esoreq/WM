clear; close all;clc;
root = '/Users/eyalsoreq/Dropbox/PhD/CURRENT/papers/WM/wm-decoding-project/results';
inputDir = [root,'/Data/MinedData'];
outDir = [root,'/report/GlobalEffects'];
mkdir(outDir)
load([root,'/Data/MinedData', '/DG_gpp.mat']);
load([root,'/report/HCL', '/DGT.mat']);
setenv('PATH', [getenv('PATH') ':/usr/local/bin'] )

tmap = [188,209,244;73,139,234;242,97,43;224,135,184;92,186,71;234,61,106]./255;
gmap = tmap([4,5,6,1,2],:);

%% mean data by clusters and processing stages
T.AAL = data.task1.BA.Properties.VariableNames(8:end)';
clusters = {'FPT','SMT','IPT','AVS','PVS'};
for ii=1:numel(clusters)
X.ba1.(clusters{ii}) = grpstats(data.task1.BA,{'subj','state'},{'mean'},'DataVars',T.AAL(T.cid==ii));
X.ba2.(clusters{ii}) = grpstats(data.task2.BA,{'subj','state'},{'mean'},'DataVars',T.AAL(T.cid==ii));
end
fn = fieldnames(X.ba1);
xx1 = zeros(19,5,2);
xx2 = zeros(16,5,2);
outDelta = table();
outRest = table();
for ii=1:numel(fn)
    f = fn{ii};
    E1 = mean(X.ba1.(f){c3nl.strDetect(X.ba1.(f).state,'Encode'),5:end},2);
    E2 = mean(X.ba2.(f){c3nl.strDetect(X.ba2.(f).state,'Encode'),5:end},2);
    M1 = mean(X.ba1.(f){c3nl.strDetect(X.ba1.(f).state,'Maintain'),5:end},2);
    M2 = mean(X.ba2.(f){c3nl.strDetect(X.ba2.(f).state,'Maintain'),5:end},2);
    [h,p,ci,stats]=ttest(E1-M1);
    outDelta = [outDelta;table({'Study1'},{f},mean(E1-M1),std(E1-M1),stats.tstat,stats.df,p,ci','VariableNames',...
    {'Study','Cluster','Mean','Std','t','df','Sig','ci'})];
    [h,p,ci,stats]=ttest(E2-M2);
    outDelta = [outDelta;table({'Study2'},{f},mean(E2-M2),std(E2-M2),stats.tstat,stats.df,p,ci','VariableNames',...
    {'Study','Cluster','Mean','Std','t','df','Sig','ci'})];
    [h,p,ci,stats]=ttest(E1);
    outRest = [outRest;table({'Study1'},{'Encode'},{f},mean(E1),std(E1),stats.tstat,stats.df,p,ci','VariableNames',...
    {'Study','Stage','Cluster','Mean','Std','t','df','Sig','ci'})];
    [h,p,ci,stats]=ttest(E2);
    outRest = [outRest;table({'Study2'},{'Encode'},{f},mean(E2),std(E2),stats.tstat,stats.df,p,ci','VariableNames',...
    {'Study','Stage','Cluster','Mean','Std','t','df','Sig','ci'})];
    [h,p,ci,stats]=ttest(M1);
    outRest = [outRest;table({'Study1'},{'Maintain'},{f},mean(M1),std(M1),stats.tstat,stats.df,p,ci','VariableNames',...
    {'Study','Stage','Cluster','Mean','Std','t','df','Sig','ci'})];
    [h,p,ci,stats]=ttest(M2);
    outRest = [outRest;table({'Study2'},{'Maintain'},{f},mean(M2),std(M2),stats.tstat,stats.df,p,ci','VariableNames',...
    {'Study','Stage','Cluster','Mean','Std','t','df','Sig','ci'})];
    xx1(:,ii,1)=E1;
    xx1(:,ii,2)=M1;
    xx2(:,ii,1)=E2;
    xx2(:,ii,2)=M2;
end

outDelta.P_FDR = get.fdr(outDelta.Sig);
outRest.P_FDR = get.fdr(outRest.Sig);
outRest = outRest(:,[1:8,10,9]);
outDelta = outDelta(:,[1:7,9,8]);
caption = ['Paired measure t-test was perfomed on the difference between activation in Encode and Delay over the FPT, SMT, IPT, AVS and PVS clusters.'];
plot.mdl2Table(outDelta,[outDir,'/GlobalDeltaEMBA.tex'],'latex',caption,'GDEMBA');

caption = ['1-sample t-test was perfomed on the mean activity of both Encode and Delay over the FPT, SMT, IPT, AVS and PVS clusters.'];
plot.mdl2Table(outRest,[outDir,'/GlobalEMBA.tex'],'latex',caption,'GEMvsRestBA');

plot.violin(xx1(:,[1,4,5],:),{{'FPT','AVS','PVS'},{'Encode','Maintain'}},hsv(12),2,[-1,3],1,'\beta')
save.pdf([outDir filesep 'violin_BA1.ai'],6,6)
plot.violin(xx2(:,[1,4,5],:),{{'FPT','AVS','PVS'},{'Encode','Maintain'}},hsv(12),2,[-1,3],0)
save.pdf([outDir filesep 'violin_BA2.ai'],6,6)

%% Now FC
clear X
X.ppi.Study1 = grpstats(data.task1.PPI,{'subj','state'},{'mean'},'DataVars',data.task1.PPI.Properties.VariableNames(7:end));%
X.ppi.Study2 = grpstats(data.task2.PPI,{'subj','state'},{'mean'},'DataVars',data.task2.PPI.Properties.VariableNames(8:end));%
tmap = [188,209,244;73,139,234;242,97,43;224,135,184;92,186,71;234,61,106]./255;
gmap = tmap([4,5,6,1,2],:);

comp = [[1:5;1:5]';nchoosek(1:5,2)];
clusters = categorical({'FPT','SMT','IPT','AVS','PVS'});

FClabels= [clusters(comp(:,1)).*clusters(comp(:,2))];
    
fn = fieldnames(X.ppi);
xx1 = zeros(19,15,2);
xx2 = zeros(16,15,2);
outDelta = table();
outRest = table();

for jj=1:numel(fn)
    f = fn{jj};
    ix = c3nl.strDetect(X.ppi.(f).state,'Encode|Maintain');
    gp = categorical(X.ppi.(f).subj(ix)).*categorical(X.ppi.(f).state(ix));
    [MU,STD,LEVELS,tmap,AB,ix_2] = get.meanconn(X.ppi.(f){ix,4:end},gp,T.number(end),T,gmap,1);
    for ii=1:size(comp,1)
        E = MU(ii,1:2:end);
        M = MU(ii,2:2:end);
        [h,p,ci,stats]=ttest(E-M);
        outDelta = [outDelta;table({f},{char(FClabels(ii))},mean(E-M),std(E-M),stats.tstat,stats.df,p,ci,'VariableNames',...
        {'Study','Cluster','Mean','Std','t','df','Sig','ci'})];
        [h,p,ci,stats]=ttest(E);
        outRest = [outRest;table({f},{'Encode'},{char(FClabels(ii))},mean(E),std(E),stats.tstat,stats.df,p,ci,'VariableNames',...
        {'Study','Stage','Cluster','Mean','Std','t','df','Sig','ci'})];
        [h,p,ci,stats]=ttest(M);
        outRest = [outRest;table({f},{'Delay'},{char(FClabels(ii))},mean(M),std(M),stats.tstat,stats.df,p,ci,'VariableNames',...
        {'Study','Stage','Cluster','Mean','Std','t','df','Sig','ci'})];
        if jj==1
            xx1(:,ii,1) = E';
            xx1(:,ii,2) = M';
        else
            xx2(:,ii,1) = E';
            xx2(:,ii,2) = M';
        end
    end
end


caption = ['Paired measure t-test was perfomed on the difference between Encode and Delay over the intra-FC clusters.'];
plot.mdl2Table(outDelta([1:5,16:20],:),[outDir,'/GlobalDeltaintraEMFC.tex'],'latex',caption,'GDEMFC');

caption = ['Paired measure t-test was perfomed on the difference between Encode and Delay over the inter-FC clusters.'];
plot.mdl2Table(outDelta([6:15,21:end],:),[outDir,'/GlobalDeltainterEMFC.tex'],'latex',caption,'GDEMFC');

caption = ['1-sample t-test was perfomed on the mean intra-connectivity of both Encode and Delay over the intra-FC clusters.'];
plot.mdl2Table(outRest([1:10,31:40],:),[outDir,'/GlobalintraEMFC.tex'],'latex',caption,'GEMvsRestFC');

caption = ['1-sample t-test was perfomed on the mean intra-connectivity of both Encode and Delay over the inter-FC clusters.'];
plot.mdl2Table(outRest([11:30,41:end],:),[outDir,'/GlobalinterEMFC.tex'],'latex',caption,'GEMvsRestFC');
c3nl.copy(c3nl.select('name','*.tex'),'/Users/eyalsoreq/Dropbox/PhD/CURRENT/papers/WM/wm-decoding-project/supplamentryMaterial')

%%
plot.violin(xx1(:,[1,4,5],:),{{'FPT','AVS','PVS'},{'Encode','Maintain'}},hsv(12),2,[-0.015,0.25],1,'\beta')
save.pdf([outDir filesep 'violin_FCintra1.ai'],6,6)
plot.violin(xx2(:,[1,4,5],:),{{'FPT','AVS','PVS'},{'Encode','Maintain'}},hsv(12),2,[-0.015,0.25],0)
save.pdf([outDir filesep 'violin_FCintra2.ai'],6,6)

plot.violin(xx1(:,[8,9,15],:),{{{'FPT','AVS'},{'FPT','PVS'},{'AVS','PVS'}},{'Encode','Maintain'}},hsv(12),2,[-0.03,0.25],1,'\beta')
save.pdf([outDir filesep 'violin_FCinter1.ai'],6,6)
plot.violin(xx2(:,[8,9,15],:),{{{'FPT','AVS'},{'FPT','PVS'},{'AVS','PVS'}},{'Encode','Maintain'}},hsv(12),2,[-0.03,0.25],0)
save.pdf([outDir filesep 'violin_FCinter2.ai'],6,6)

%% motor effects
% 
% ix = ~c3nl.strDetect(X.ppi.Study1.state,'Cue');
% gp = categorical(X.ppi.Study1.subj(ix)).*categorical(X.ppi.Study1.state(ix));
% [MU,STD,LEVELS,tmap,AB,ix2] = get.meanconn(X.ppi.Study1{ix,4:end},gp,T.number(end),T,gmap,1);
% S1 = MU(2,:)'.*dummyvar(categorical(X.ppi.Study1.state(ix)));
% 
% ix = ~c3nl.strDetect(X.ppi.Study2.state,'Cue');
% gp = categorical(X.ppi.Study2.subj(ix)).*categorical(X.ppi.Study2.state(ix));
% [MU,STD,LEVELS,tmap,AB,ix2] = get.meanconn(X.ppi.Study2{ix,4:end},gp,T.number(end),T,gmap,1);
% S2 = MU(2,:)'.*dummyvar(categorical(X.ppi.Study2.state(ix)));
% 
% xx1p = zeros(nnz(S1(:,2)),1,3);
% for ii=1:3;xx1p(:,1,ii)=S1(S1(:,ii)~=0,ii);end
% plot.violin(xx1p,{{'SMT'},{'Encode','Maintain','probe'}},hsv(12),2,[-.1,.4],0)
% save.pdf([outDir filesep 'violin_SMTFC1.ai'],3,6)
% 
% xx2p = zeros(nnz(S2(:,2)),1,3);
% for ii=1:3;xx2p(:,1,ii)=S2(S2(:,ii)~=0,ii);end
% plot.violin(xx2p,{{'SMT'},{'Encode','Maintain','probe'}},hsv(12),2,[-.1,.4],0)
% save.pdf([outDir filesep 'violin_SMTFC2.ai'],3,6)
% 
% 
% X.ba1 = grpstats(data.task1.BA,{'subj','state'},{'mean'},'DataVars',T.name(T.cid==2));% SMT
% X.ba2 = grpstats(data.task2.BA,{'subj','state'},{'mean'},'DataVars',T.name(T.cid==2));% SMT
% ix = ~c3nl.strDetect(cellstr(X.ba1.state),'Cue');
% S1 = mean(X.ba1{ix,5:end},2).*dummyvar(categorical(X.ba1.state(ix)));
% xx1 = zeros(nnz(S1(:,2)),1,3);
% for ii=1:3;xx1(:,1,ii)=S1(S1(:,ii)~=0,ii);end
% plot.violin(xx1,{{'SMT'},{'Encode','Maintain','probe'}},hsv(12),2,[-1,3.5],0)
% save.pdf([outDir filesep 'violin_SMTBA1.ai'],3,6)
% 
% ix = ~c3nl.strDetect(cellstr(X.ba2.state),'Cue');
% S2 = mean(X.ba2{ix,5:end},2).*dummyvar(categorical(X.ba2.state(ix)));
% xx2 = zeros(nnz(S2(:,2)),1,3);
% for ii=1:3;xx2(:,1,ii)=S2(S2(:,ii)~=0,ii);end
% plot.violin(xx2,{{'SMT'},{'Encode','Maintain','probe'}},hsv(12),2,[-1,3.5],0)
% save.pdf([outDir filesep 'violin_SMTBA2.ai'],3,6)
% 
% ix1 = ~c3nl.strDetect(cellstr(X.ba1.state),'Cue');
% ix2 = ~c3nl.strDetect(cellstr(X.ba2.state),'Cue');
% 
% gp1 = categorical(X.ppi.Study1.subj(ix1)).*categorical(X.ppi.Study1.state(ix1));
% MU1 =  get.meanconn(X.ppi.Study1{ix1,4:end},gp1,T.number(end),T,gmap,1);
% gp2 = categorical(X.ppi.Study2.subj(ix2)).*categorical(X.ppi.Study2.state(ix2));
% MU2 =  get.meanconn(X.ppi.Study2{ix2,4:end},gp2,T.number(end),T,gmap,1);
% 
% 
% LongY = [table(categorical([repmat({'study1'},nnz(ix1)*2,1);repmat({'study2'},nnz(ix2)*2,1)]),'VariableNames',{'Study'}),...
%     table(categorical([cellstr(X.ba1.subj(ix1,:));X.ppi.Study1.subj(ix1);cellstr(X.ba2.subj(ix2,:));X.ppi.Study2.subj(ix2)]),'VariableNames',{'Subj'}),...
%     table(categorical([cellstr(X.ba1.state(ix1,:));X.ppi.Study1.state(ix1);cellstr(X.ba2.state(ix2,:));X.ppi.Study2.state(ix2)]),'VariableNames',{'Stage'}),...
%     table(categorical([repmat({'BA'},nnz(ix1),1);repmat({'FC'},nnz(ix1),1);repmat({'BA'},nnz(ix2),1);repmat({'FC'},nnz(ix2),1)]),'VariableNames',{'Metric'}),...
%     array2table([c3nl.scale(mean(X.ba1{ix1,5:end},2));c3nl.scale(MU1(2,:))';c3nl.scale(mean(X.ba2{ix2,5:end},2));c3nl.scale(MU2(2,:))'],'VariableNames',{'Ys'}),...
%     array2table([zscore(mean(X.ba1{ix1,5:end},2));zscore(MU1(2,:))';zscore(mean(X.ba2{ix2,5:end},2));zscore(MU2(2,:))'],'VariableNames',{'Yz'}),...
%     array2table([mean(X.ba1{ix1,5:end},2);MU1(2,:)';mean(X.ba2{ix2,5:end},2);MU2(2,:)'],'VariableNames',{'Y'})];
% 
% lme = fitlme(LongY, 'Y ~ 1 + Study+Stage+Metric+(1|Subj)','DummyVarCoding','effects');
% plot.mdl2Table(lme.Coefficients,[outDir filesep 'SMT_lmestats.tex'],'latex','SMT increases as a function of probe stage regardless of metric across studies','SMTlme')
% plot.mdl2Table(anova(lme),[outDir filesep 'SMT_anovastats.tex'],'latex','SMT probe stage ANOVA(lme) main effects ','SMTlmeAnova')
% 


%% load effects take two
% BOLD activity Load data
ix1 = c3nl.strDetect(data.task1.BA.state,'Maintain');
ix2 = c3nl.strDetect(data.task2.BA.state,'Maintain');

T.AAL = data.task1.BA.Properties.VariableNames(8:end)';
cl = struct('FPT',1,'AVS',4,'PVS',5);
fn = fieldnames(cl);
X.ba1 = table();
X.ba2 = table();
for ii = 1:numel(fn)
    tmp = grpstats(data.task1.BA(ix1,:),{'subj','load'},{'mean'},'DataVars',T.AAL(ismember(T.cid,cl.(fn{ii}))));% frontoparietal
    X.ba1 = [X.ba1;table(tmp.subj,tmp.load,repmat(fn{ii},height(tmp),1),repmat('BA',height(tmp),1),mean(tmp{:,4:end},2),'VariableNames',{'subj','load','cluster','metric','Y'})];
    tmp = grpstats(data.task2.BA(ix2,:),{'subj','load'},{'mean'},'DataVars',T.AAL(ismember(T.cid,cl.(fn{ii}))));% frontoparietal
    X.ba2 = [X.ba2;table(tmp.subj,tmp.load,repmat(fn{ii},height(tmp),1),repmat('BA',height(tmp),1),mean(tmp{:,4:end},2),'VariableNames',{'subj','load','cluster','metric','Y'})];
end
X.ba1.load = categorical(X.ba1.load,{'Low','Med','High'},{'Low','Med','High'});
X.ba2.load = categorical(X.ba2.load,{'Low','High'},{'Low','High'});
ix1 = c3nl.strDetect(data.task1.PPI.state,'Maintain');
ix2 = c3nl.strDetect(data.task2.PPI.state,'Maintain');
X.ppi1 = grpstats(data.task1.PPI(ix1,:),{'subj','load'},{'mean'},'DataVars',data.task1.PPI.Properties.VariableNames(7:end));%
X.ppi2 = grpstats(data.task2.PPI(ix2,:),{'subj','load'},{'mean'},'DataVars',data.task2.PPI.Properties.VariableNames(8:end));%

gp = categorical(X.ppi1.subj.*categorical(X.ppi1.load,{'Low','Med','High'},{'Low','Med','High'}));
MU1 = get.meanconn(X.ppi1{:,4:end},gp,T.number(end),T,gmap,1);
gp = categorical(X.ppi2.subj.*categorical(X.ppi2.load,{'Low','High'},{'Low','High'}));
MU2 = get.meanconn(X.ppi2{:,4:end},gp,T.number(end),T,gmap,1);
pairs = {'FP2AV';'FP2PV';'AV2PV'};

X.inter1 = table(categorical([X.ppi1.subj;X.ppi1.subj;X.ppi1.subj]),...
                  categorical([X.ppi1.load;X.ppi1.load;X.ppi1.load],{'Low','Med','High'},{'Low','Med','High'}),...
                  repelem(pairs,size(MU1,2),1),repmat('Inter',size(MU1,2)*3,1),...
                  [MU1(8,:),MU1(9,:),MU1(15,:)]','VariableNames',{'subj','load','cluster','metric','Y'});

X.inter2 = table(categorical([X.ppi2.subj;X.ppi2.subj;X.ppi2.subj]),...
                  categorical([X.ppi2.load;X.ppi2.load;X.ppi2.load],{'Low','Med','High'},{'Low','Med','High'}),...
                  repelem(pairs,size(MU2,2),1),repmat('Inter',size(MU2,2)*3,1),...
                  [MU2(8,:),MU2(9,:),MU2(15,:)]','VariableNames',{'subj','load','cluster','metric','Y'});

X.intra1 = table(categorical([X.ppi1.subj;X.ppi1.subj;X.ppi1.subj]),...
                  categorical([X.ppi1.load;X.ppi1.load;X.ppi1.load],{'Low','Med','High'},{'Low','Med','High'}),...
                  repelem(fn,size(MU1,2),1),repmat('Intra',size(MU1,2)*3,1),...
                  [MU1(1,:),MU1(4,:),MU1(5,:)]','VariableNames',{'subj','load','cluster','metric','Y'});

X.intra2 =table(categorical([X.ppi2.subj;X.ppi2.subj;X.ppi2.subj]),...
                  categorical([X.ppi2.load;X.ppi2.load;X.ppi2.load],{'Low','High'},{'Low','High'}),...
                  repelem(fn,size(MU2,2),1),repmat({'Intra'},size(MU2,2)*3,1),...
                  [MU2(1,:),MU2(4,:),MU2(5,:)]','VariableNames',{'subj','load','cluster','metric','Y'});

% standerdize each group indpendently 

X.ba1.YZ = zscore(X.ba1.Y);
X.ba2.YZ = zscore(X.ba2.Y);
X.inter1.YZ = zscore(X.inter1.Y);
X.inter2.YZ = zscore(X.inter2.Y);
X.intra1.YZ = zscore(X.intra1.Y);
X.intra2.YZ = zscore(X.intra2.Y);
fn = {'ba','inter','intra'};
study1 = table();
study2 = table();
for ii=1:numel(fn)
tmp = grpstats(X.([fn{ii} '1']),{'subj','load','metric'},{'mean'},'DataVars','YZ');%
study1 = [study1;table(categorical(cellstr(tmp.subj)),tmp.load,categorical(cellstr(tmp.metric)),tmp.mean_YZ,'VariableNames',{'subj','load','metric','Y'})];
tmp = grpstats(X.([fn{ii} '2']),{'subj','load','metric'},{'mean'},'DataVars','YZ');%
study2 = [study2;table(categorical(cellstr(tmp.subj)),tmp.load,categorical(cellstr(tmp.metric)),tmp.mean_YZ,'VariableNames',{'subj','load','metric','Y'})];
end

gp = study1.metric.*study1.load;
wide1 = unstack(table(study1.subj,gp,study1.Y,'VariableNames',{'subj','gp','Y'}),'Y','gp');
factors = convert.cell2Table(cellfun(@(x) strsplit(x,' '),categories(gp)','un',0)');
factors.Properties.VariableNames = {'metric','load'};
mes = wide1.Properties.VariableNames(2:end);
rm = fitrm(wide1,sprintf('%s-%s~1',mes{1},mes{end}),'WithinDesign',factors);
mauchly(rm)
epsilon(rm)
tbl_RMABA1 = ranova(rm,'WithinModel',strjoin({'metric','load'},'*'));

rm1 = fitrm(wide1(:,1:4),sprintf('%s-%s~1',mes{1},mes{3}),'WithinDesign',factors(1:3,2));
mauchly(rm1)
epsilon(rm1)
tbl_RMABA1oneway = ranova(rm1,'WithinModel','load');
rm2 = fitrm(wide1(:,[1,5:7]),sprintf('%s-%s~1',mes{4},mes{6}),'WithinDesign',factors(1:3,2));
mauchly(rm2)
epsilon(rm2)
tbl_RMAinter1oneway = ranova(rm2,'WithinModel','load');
tbl_RMAinter1oneway_PH = c3nl.addt(multcompare(rm2,'load','ComparisonType' ,'bonferroni'));
rm3 = fitrm(wide1(:,[1,8:10]),sprintf('%s-%s~1',mes{7},mes{9}),'WithinDesign',factors(1:3,2));
mauchly(rm3)
epsilon(rm3)
tbl_RMAintra1oneway = ranova(rm3,'WithinModel','load');
tbl_RMAintra1oneway_PH = c3nl.addt(multcompare(rm3,'load','ComparisonType' ,'bonferroni'));

caption = ['Repeated measures two way ANOVA was used to interogate the effects of load over metric in study one.'];
plot.mdl2Table(c3nl.roundTable(tbl_RMABA1(:,[1:4,6]),4),[outDir filesep 'tbl_RMABA1.tex'],'latex',caption,'tbl_RMABA1');
caption = ['Repeated measures one way ANOVA was used to interogate the effects of load for activation in study one.'];
plot.mdl2Table(c3nl.roundTable(tbl_RMABA1oneway(:,[1:4,6]),4),[outDir filesep 'tbl_RMABA1oneway.tex'],'latex',caption,'tbl_RMABA1');
caption = ['Repeated measures one way ANOVA was used to interogate the effects of load for inter-connectivity in study one.'];
plot.mdl2Table(c3nl.roundTable(tbl_RMAinter1oneway(:,[1:5]),4),[outDir filesep 'tbl_RMAinter1oneway.tex'],'latex',caption,'tbl_RMABA1');
caption = ['Repeated measures one way ANOVA was used to interogate the effects of load for intra-connectivity in study one.'];
plot.mdl2Table(c3nl.roundTable(tbl_RMAintra1oneway(:,[1:5]),4),[outDir filesep 'tbl_RMAintra1oneway.tex'],'latex',caption,'tbl_RMABA1');
caption = ['Multicomparison post hoc test (bonferroni corrected) was used to interogate the significant main effect of load found in the Repeated measures ANOVA for inter-connectivity.'];
plot.mdl2Table(c3nl.roundTable(tbl_RMAinter1oneway_PH,4),[outDir filesep 'tbl_RMAinter1oneway_PH.tex'],'latex',caption,'tbl_RMABA1_PH_load');
caption = ['Multicomparison post hoc test (bonferroni corrected) was used to interogate the significant main effect of load found in the Repeated measures ANOVA for intra-connectivity.'];
plot.mdl2Table(c3nl.roundTable(tbl_RMAintra1oneway_PH,4),[outDir filesep 'tbl_RMAintra1oneway_PH.tex'],'latex',caption,'tbl_RMABA1_PH_load');


gp = study2.metric.*study2.load;
wide2 = unstack(table(study2.subj,gp,study2.Y,'VariableNames',{'subj','gp','Y'}),'Y','gp');
factors = convert.cell2Table(cellfun(@(x) strsplit(x,' '),cellstr(unique(gp))','un',0)');
factors.Properties.VariableNames = {'metric','load'};
mes = wide2.Properties.VariableNames(2:end);
rm = fitrm(wide2,sprintf('%s-%s~1',mes{1},mes{end}),'WithinDesign',factors);
mauchly(rm)
epsilon(rm)
tbl_RMABA2 = ranova(rm,'WithinModel',strjoin({'metric','load'},'*'));




tbl_RMABA2_PH_load = c3nl.addt(multcompare(rm,'load','ComparisonType' ,'bonferroni'));
tbl_RMABA2_PH_inter = c3nl.addt(multcompare(rm,'load','By', 'metric','ComparisonType' ,'bonferroni'));

caption = ['Repeated measures ANOVA was used to interogate the effects of load over stage for activation in study one.'];
plot.mdl2Table(c3nl.roundTable(tbl_RMABA2(:,[1:4,6]),4),[outDir filesep 'tbl_RMABA2.tex'],'latex',caption,'tbl_RMABA2');
caption = ['Multicomparison post hoc test (bonferroni corrected) was used to interogate the significant main effect of load found in the Repeated measures ANOVA.'];
plot.mdl2Table(c3nl.roundTable(tbl_RMABA2_PH_load,4),[outDir filesep 'tbl_RMABA2_PH_load.tex'],'latex',caption,'tbl_RMABA2_PH_load');
caption = ['Multicomparison post hoc test (bonferroni corrected) was used to interogate the significant interaction of load by stage found in the Repeated measures ANOVA.'];
plot.mdl2Table(c3nl.roundTable(tbl_RMABA2_PH_inter,4),[outDir filesep 'tbl_RMABA2_PH_inter.tex'],'latex',caption,'tbl_RMABA2_PH_inter');


c3nl.copy(c3nl.select('name','*.tex'),'/Users/eyalsoreq/Dropbox/PhD/CURRENT/papers/WM/wm-decoding-project/supplamentryMaterial')
%% plot figures 


clear X

T.AAL = data.task1.BA.Properties.VariableNames(8:end)';
X.ba1 = grpstats(data.task1.BA,{'subj','state','load'},{'mean'},'DataVars',T.AAL(ismember(T.cid,[1,4,5])));
X.ba1.load = categorical(X.ba1.load,{'Low','Med','High'},{'Low','Med','High'});
tmp = unstack(table(X.ba1.subj,X.ba1.state,X.ba1.load,mean(X.ba1{:,5:end},2),'VariableNames',{'subj','state','load','Y'}),'Y','load');
X.ba1 = grpstats(tmp,'state',{'mean','std'},'DataVars',{'Low','Med','High'});
X.ba2 = grpstats(data.task2.BA,{'subj','state','load'},{'mean'},'DataVars',T.AAL(ismember(T.cid,[1,4,5])));
X.ba2.load = categorical(X.ba2.load,{'Low','High'},{'Low','High'});
tmp = unstack(table(X.ba2.subj,X.ba2.state,X.ba2.load,mean(X.ba2{:,5:end},2),'VariableNames',{'subj','state','load','Y'}),'Y','load');
X.ba2 = grpstats(tmp,'state',{'mean','std'},'DataVars',{'Low','High'});

tmp = grpstats(data.task1.PPI,{'subj','state','load'},{'mean'},'DataVars',data.task1.PPI.Properties.VariableNames(7:end));%
tmp.load = categorical(tmp.load,{'Low','Med','High'},{'Low','Med','High'});
gp = categorical(tmp.subj).*tmp.state.*tmp.load;
MU1 = get.meanconn(tmp{:,5:end},gp,T.number(end),T,gmap,1);
gp = convert.cell2Table(arrayfun(@(x) strsplit(char(x),' '), gp,'un',0),{'subj','state','load'});
tmp = unstack([gp,table(mean(MU1([1,4:5],:)',2),'VariableNames',{'Y'})],'Y','load');
X.ppi1_intra = grpstats(tmp,'state',{'mean','std'},'DataVars',{'Low','Med','High'});
tmp = unstack([gp,table(mean(MU1([8:9,15],:)',2),'VariableNames',{'Y'})],'Y','load');
X.ppi1_inter = grpstats(tmp,'state',{'mean','std'},'DataVars',{'Low','Med','High'});


tmp = grpstats(data.task2.PPI,{'subj','state','load'},{'mean'},'DataVars',data.task2.PPI.Properties.VariableNames(8:end));%
tmp.load = categorical(tmp.load,{'Low','High'},{'Low','High'});
gp = categorical(tmp.subj).*tmp.state.*tmp.load;
MU2 = get.meanconn(tmp{:,5:end},gp,T.number(end),T,gmap,1);
gp = convert.cell2Table(arrayfun(@(x) strsplit(char(x),' '), gp,'un',0),{'subj','state','load'});
tmp = unstack([gp,table(mean(MU2([1,4:5],:)',2),'VariableNames',{'Y'})],'Y','load');
X.ppi2_intra = grpstats(tmp,'state',{'mean','std'},'DataVars',{'Low','High'});
tmp = unstack([gp,table(mean(MU2([8:9,15],:)',2),'VariableNames',{'Y'})],'Y','load');
X.ppi2_inter = grpstats(tmp,'state',{'mean','std'},'DataVars',{'Low','High'});



cmap =[0,176,254;0,104,253;0,40,168]./256;

fn = fieldnames(X);
figure(1);ax=gca;
for ii=1:numel(fn)
    cla(ax)
    if c3nl.strDetect(fn(ii),'1')
        plot.splines2d(1:4,X.(fn{ii}){1:4,3:2:7}',ax,cmap,X.(fn{ii}).state(1:4),X.(fn{ii}){1:4,4:2:8}'./sqrt(X.(fn{ii}).GroupCount(1)),4,[],0.175)
    else 
        plot.splines2d(1:4,X.(fn{ii}){1:4,3:2:5}',ax,cmap([1,3],:),X.(fn{ii}).state(1:4),X.(fn{ii}){1:4,4:2:6}'./sqrt(X.(fn{ii}).GroupCount(1)),4,[],0.175)
    end
    save.pdf([outDir filesep 'spline_' fn{ii} '.ai'],5,2)
    
end

