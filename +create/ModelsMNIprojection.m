clear; close all;clc;
root = '/Users/eyalsoreq/Dropbox/PhD/CURRENT/papers/WM/wm-decoding-project/results/report';
outDir = [root,'/SparseProjection'];
mkdir(outDir);

%% start with DG for the paper 

load /Users/eyalsoreq/Dropbox/PhD/CURRENT/papers/WM/wm-decoding-project/results/Data/MinedData/DG_gpp.mat  
Xt_ds = grpstats(data.task1.PPI,{'subj','run','domain','state'},{'mean'},'DataVars',data.task1.PPI.Properties.VariableNames(7:end));% aggregating over dfact
XT_ds = grpstats(data.task2.PPI,{'subj','run','domain','state'},{'mean'},'DataVars',data.task2.PPI.Properties.VariableNames(8:end));% aggregating over dfact
ixT = ~c3nl.strDetect(XT_ds.state,'Cue');ixt = ~c3nl.strDetect(Xt_ds.state,'Cue');
xTm = grpstats(XT_ds(ixT,:),{'subj','run','domain'},{'mean'},'DataVars',XT_ds.Properties.VariableNames(6:end));% aggregating over state
Xt = Xt_ds{ixt,6:end};% training 
XT = XT_ds{ixT,6:end};% test with stage
xxT = xTm{:,5:end};% test with stage avreged 
Yt = categorical(Xt_ds.domain(ixt));YT = categorical(XT_ds.domain(ixT),categories(Yt));yT = categorical(xTm.domain,categories(Yt));
St = categorical(Xt_ds.state(ixt));ST = categorical(XT_ds.state(ixT));
addpath /Users/eyalsoreq/Dropbox/PhD/CURRENT/code/c3nl_fusion/3rdParty/glmnet_matlab      
op = glmnetSet;
op.mtype = 'ungrouped';
op.ltype= 'Newton';
f = @(t,mx,mn,a) exp((t./numel(t)).*log(a*mn/mx))./mn;
op.lambda = f(1:100,size(Xt,2),numel(unique(Yt)),1);
rng(22041975)
cvmdl = ml.fit.glmnet(Xt,Yt,'nfolds',5,'opt',op);
mdl = cvmdl.glmnet_fit;
mdl.cvmdl = cvmdl;
[~,yp]=ml.predict.glmnet(mdl,Xt);
idx=plot.scree(squeeze(mean(min(ml.score.crossEntropy(yp),[],2))),1);
mdl.fs = [mdl.beta{1}(:,idx),mdl.beta{2}(:,idx),mdl.beta{3}(:,idx)];
mdl.idx = idx;
df = nnz(mdl.fs);
fsx = any(mdl.fs,2);
name = 'DG';
cmap = [151,84,204;129,198,121;243,184,99]./256;% frac/Num/Pos
gmap = repelem(cmap,3,1).^repmat([0.9;2;5],3,1);% frac/Num/Pos
[yc,yp] =ml.predict.glmnet(mdl,Xt,[],mdl.lambda(idx));
plot.radviz([Yt.*St,yc.*St],yp,categories(Yt),cmap,gmap,4);
save.pdf(sprintf('%s%sfs_%s_RVtrain.ai',outDir,filesep,name),12,12);
plot.confusiongrad(Yt,yc,cellstr(unique(Yt)),sprintf('df=%i',df),600,600,0.05,2,[0 0 0 1]);
save.pdf(sprintf('%s%sfs_%s_CMtrain.ai',outDir,filesep,name),6,6);

[yc,yp] =ml.predict.glmnet(mdl,XT,[],mdl.lambda(idx));
plot.radviz([YT.*ST,yc.*ST],yp,categories(Yt),cmap,gmap,4);
save.pdf(sprintf('%s%sfs_%s_RVtest.ai',outDir,filesep,name),12,12);
plot.confusiongrad(YT,yc,cellstr(unique(Yt)),sprintf('df=%i',df),600,600,0.05,2,[0 0 0 1]);
save.pdf(sprintf('%s%sfs_%s_CMtest.ai',outDir,filesep,name),6,6);

[yc,yp] =ml.predict.glmnet(mdl,xxT,[],mdl.lambda(idx));
plot.radviz([yT,yc],yp,categories(Yt),cmap,cmap,4);
save.pdf(sprintf('%s%sfs_%s_RVmutest.ai',outDir,filesep,name),12,12);
plot.confusiongrad(yT,yc,cellstr(unique(Yt)),sprintf('df=%i',df),600,600,0.05,2,[0 0 0 1]);
save.pdf(sprintf('%s%sfs_%s_CMmutest.ai',outDir,filesep,name),6,6);


[~,G]=load.vol('/Users/eyalsoreq/Dropbox/PhD/CURRENT/papers/WM/wm-decoding-project/results/Data/conjunction/task1/spmT_0006_Stages_cluster_corrected.nii');
cmap = [151,84,204;129,198,121;243,184,99]./256;% frac/Num/Pos

tmap = hsv(360);
fmap = get.cmap([tmap(220,:);tmap(220,:).^.5;.85,.85,.85;cmap(1,:);cmap(1,:).^2]);
nmap = get.cmap([tmap(220,:);tmap(220,:).^.5;.85,.85,.85;cmap(2,:);cmap(2,:).^3]);
pmap = get.cmap([tmap(220,:);tmap(220,:).^.5;.85,.85,.85;cmap(3,:);cmap(3,:).^3]);

TT = atlas.toTable(data.roi.roi.L,G);

plot.projectedAdj(mdl.fs,TT,struct('or','a','cmap',cmap));
save.pdf(sprintf('%s%sfs_%s_mdl_Axial.ai',outDir,filesep,name),6,6);
plot.projectedAdj(mdl.fs,TT,struct('or','rs','cmap',cmap));
save.pdf(sprintf('%s%sfs_%s_mdl_Saggital.ai',outDir,filesep,name),6,6);
plot.projectedAdj(mdl.fs,TT,struct('or','c','cmap',cmap));
save.pdf(sprintf('%s%sfs_%s_mdl_Coronal.ai',outDir,filesep,name),6,6);

plot.multiCAdj(TT,mdl.fs,cmap,10,2);
save.pdf(sprintf('%s%sfs_%s_Circ.ai',outDir,filesep,name),15,15);

cmap = [224,135,184;92,186,71;234,61,106;188,209,244;73,139,234;242,97,43]./255;

load('/Users/eyalsoreq/Dropbox/PhD/CURRENT/papers/WM/wm-decoding-project/results/report/HCL/DGT.mat')
domains = {'Frac','Num','Pos'};
for ii=1:3
    plot.vane(mdl.fs(:,ii),T,cmap,1,0,600,600);
    save.pdf(sprintf('%s%sfs_%s_vane_%s.ai',outDir,filesep,name,domains{ii}),6,6);
end
outT = table();
[ix,m] =find(mdl.fs);
weights= mdl.fs(ix,:);
idx=find(triu(true(height(T)),1));
for ii=1:3
    [I,J]=ind2sub([height(T),height(T)],idx(ix(m==ii))); 
    w = weights(m==ii,ii);
    labels = [TT.name_AAL2(I) TT.name_AAL2(J)];
    tmp = XT(:,ix(m==ii));
    tn = tmp(YT=='Num',:);
    tp = tmp(YT~='Num',:);
    for jj=1:size(tp,2)
        [h,p,ci,stats]=ttest2(tn(:,jj),tp(:,jj));
        outT = [outT;table(labels(jj,1),labels(jj,2),domains(ii),w(jj), mean(tn(:,jj)-tp(:,jj)),std(tn(:,jj)-tp(:,jj)),stats.tstat,stats.df,p,ci','VariableNames',...
        {'node_1','node_2','domain','weight','Mean_diff','Std_diff','t','df','Sig','ci'})];
    end
end
outT.p_FDR = get.fdr(outT.Sig);
[I,J]=ind2sub([height(T),height(T)],idx(ix));
cmap = [151,84,204;129,198,121;243,184,99]./256;% frac/Num/Pos
caption = ['Two sample t-test was perfomed on functional connectivity measures extracted from the three class domains models.'];
plot.mdl2Table(c3nl.roundTable(outT(:,[1:7,9,end]),3),[outDir,'/TTestFinal.tex'],'latex',caption,'TTestFinal');

[v,ix] =max(mdl.fs);

plot.bar(XT(:,ix),YT,{'12-13','8-23','14-31'},cmap(2:3,:),1,[0,0.16],16)
save.pdf(sprintf('%s%sfs_%s_barXT.ai',outDir,filesep,name),6,6);
plot.bar(Xt(:,ix),Yt,{'12-13','8-23','14-31'},cmap,1,[0,0.16],19)
save.pdf(sprintf('%s%sfs_%s_barXt.ai',outDir,filesep,name),6,6);
c3nl.copy(c3nl.select('name','*.tex'),'/Users/eyalsoreq/Dropbox/PhD/CURRENT/papers/WM/wm-decoding-project/supplamentryMaterial')



%% do the same for SG

load /Users/eyalsoreq/Dropbox/PhD/CURRENT/papers/WM/wm-decoding-project/results/Data/MinedData/SG_gpp.mat  
Xt_ds = grpstats(data.task1.PPI,{'subj','run','domain','state'},{'mean'},'DataVars',data.task1.PPI.Properties.VariableNames(7:end));% aggregating over dfact
XT_ds = grpstats(data.task2.PPI,{'subj','run','domain','state'},{'mean'},'DataVars',data.task2.PPI.Properties.VariableNames(8:end));% aggregating over dfact
ixT = ~c3nl.strDetect(XT_ds.state,'Cue');ixt = ~c3nl.strDetect(Xt_ds.state,'Cue');
xTm = grpstats(XT_ds(ixT,:),{'subj','run','domain'},{'mean'},'DataVars',XT_ds.Properties.VariableNames(6:end));% aggregating over state
Xt = Xt_ds{ixt,6:end};% training 
XT = XT_ds{ixT,6:end};% test with stage
xxT = xTm{:,5:end};% test with stage avreged 
Yt = categorical(Xt_ds.domain(ixt));YT = categorical(XT_ds.domain(ixT),categories(Yt));yT = categorical(xTm.domain,categories(Yt));
St = categorical(Xt_ds.state(ixt));ST = categorical(XT_ds.state(ixT));
addpath /Users/eyalsoreq/Dropbox/PhD/CURRENT/code/c3nl_fusion/3rdParty/glmnet_matlab      
op = glmnetSet;
op.mtype = 'ungrouped';
op.ltype= 'Newton';
f = @(t,mx,mn,a) exp((t./numel(t)).*log(a*mn/mx))./mn;

op.lambda = f(1:100,size(Xt,2),numel(unique(Yt)),1);
rng(22041975)
cvmdl = ml.fit.glmnet(Xt,Yt,'nfolds',5,'opt',op);
mdl = cvmdl.glmnet_fit;
mdl.cvmdl = cvmdl;
[~,yp]=ml.predict.glmnet(mdl,Xt);
idx=plot.scree(squeeze(sum(min(-log(yp),[],2))),1);
mdl.fs = [mdl.beta{1}(:,idx),mdl.beta{2}(:,idx),mdl.beta{3}(:,idx)];
mdl.idx = idx;
df = nnz(mdl.fs);
fsx = any(mdl.fs,2);
name = 'SG';
[yc,yp] =ml.predict.glmnet(mdl,Xt,[],mdl.lambda(idx));
plot.radviz([Yt.*St,yc.*St],yp,categories(Yt),cmap,gmap,4);
save.pdf(sprintf('%s%sfs_%s_RVtrain.ai',outDir,filesep,name),12,12);
plot.confusiongrad(Yt,yc,cellstr(unique(Yt)),sprintf('df=%i',df),600,600,0.05,2,[0 0 0 1]);
save.pdf(sprintf('%s%sfs_%s_CMtrain.ai',outDir,filesep,name),6,6);

[yc,yp] =ml.predict.glmnet(mdl,XT,[],mdl.lambda(idx));
plot.radviz([YT.*ST,yc.*ST],yp,categories(Yt),cmap,gmap,4);
save.pdf(sprintf('%s%sfs_%s_RVtest.ai',outDir,filesep,name),12,12);
plot.confusiongrad(YT,yc,cellstr(unique(Yt)),sprintf('df=%i',df),600,600,0.05,2,[0 0 0 1]);
save.pdf(sprintf('%s%sfs_%s_CMtest.ai',outDir,filesep,name),6,6);

[yc,yp] =ml.predict.glmnet(mdl,xxT,[],mdl.lambda(idx));
plot.radviz([yT,yc],yp,categories(Yt),cmap,cmap,4);
save.pdf(sprintf('%s%sfs_%s_RVmutest.ai',outDir,filesep,name),12,12);
plot.confusiongrad(yT,yc,cellstr(unique(Yt)),sprintf('df=%i',df),600,600,0.05,2,[0 0 0 1]);
save.pdf(sprintf('%s%sfs_%s_CMmutest.ai',outDir,filesep,name),6,6);

T = atlas.toTable(data.roi.roi.L,G);
plot.projectedAdj(mdl.fs,T,struct('or','a','cmap',cmap));
save.pdf(sprintf('%s%sfs_%s_mdl_Axial.ai',outDir,filesep,name),6,6);
plot.projectedAdj(mdl.fs,T,struct('or','rs','cmap',cmap));
save.pdf(sprintf('%s%sfs_%s_mdl_Saggital.ai',outDir,filesep,name),6,6);
plot.projectedAdj(mdl.fs,T,struct('or','c','cmap',cmap));
save.pdf(sprintf('%s%sfs_%s_mdl_Coronal.ai',outDir,filesep,name),6,6);

plot.multiCAdj(T,mdl.fs,cmap,10,2);
save.pdf(sprintf('%s%sfs_%s_Circ.ai',outDir,filesep,name),6,6);

cmap = [224,135,184;92,186,71;234,61,106;188,209,244;73,139,234;242,97,43]./255;

%% do the same for WB

load ('/Users/eyalsoreq/Dropbox/PhD/CURRENT/papers/WM/wm-decoding-project/results/Data/MinedData/shen268_gpp.mat')
Xt_ds = grpstats(data.task1.PPI,{'subj','run','domain','state'},{'mean'},'DataVars',data.task1.PPI.Properties.VariableNames(7:end));% aggregating over dfact
XT_ds = grpstats(data.task2.PPI,{'subj','run','domain','state'},{'mean'},'DataVars',data.task2.PPI.Properties.VariableNames(8:end));% aggregating over dfact
ixT = ~c3nl.strDetect(XT_ds.state,'Cue');ixt = ~c3nl.strDetect(Xt_ds.state,'Cue');
xTm = grpstats(XT_ds(ixT,:),{'subj','run','domain'},{'mean'},'DataVars',XT_ds.Properties.VariableNames(6:end));% aggregating over state
Xt = Xt_ds{ixt,6:end};% training 
XT = XT_ds{ixT,6:end};% test with stage
xxT = xTm{:,5:end};% test with stage avreged 
Yt = categorical(Xt_ds.domain(ixt));YT = categorical(XT_ds.domain(ixT),categories(Yt));yT = categorical(xTm.domain,categories(Yt));
St = categorical(Xt_ds.state(ixt));ST = categorical(XT_ds.state(ixT));
op = glmnetSet;
op.mtype = 'ungrouped';
op.ltype= 'Newton';
f = @(t,mx,mn,tx) mn.*exp((t./tx).*log(mx/mn));
op.lambda = f(1:100,1/(size(Xt,2)/2),1/numel(unique(Yt)),100);
rng(22041975)
cvmdl = ml.fit.glmnet(Xt,Yt,'nfolds',5,'opt',op);
mdl = cvmdl.glmnet_fit;
mdl.cvmdl = cvmdl;
[~,yp]=ml.predict.glmnet(mdl,Xt);
idx=plot.scree(squeeze(mean(min(ml.score.crossEntropy(yp),[],2))),1);
mdl.fs = [mdl.beta{1}(:,idx),mdl.beta{2}(:,idx),mdl.beta{3}(:,idx)];
mdl.idx = idx;
df = nnz(mdl.fs);
fsx = any(mdl.fs,2);
name = 'WB';
[yc,yp] =ml.predict.glmnet(mdl,Xt,[],mdl.lambda(idx));
plot.radviz([Yt.*St,yc.*St],yp,categories(Yt),cmap,gmap,4);
save.pdf(sprintf('%s%sfs_%s_RVtrain.ai',outDir,filesep,name),12,12);
plot.confusiongrad(Yt,yc,cellstr(unique(Yt)),sprintf('df=%i',df),600,600,0.05,2,[0 0 0 1]);
save.pdf(sprintf('%s%sfs_%s_CMtrain.ai',outDir,filesep,name),6,6);

[yc,yp] =ml.predict.glmnet(mdl,XT,[],mdl.lambda(idx));
plot.radviz([YT.*ST,yc.*ST],yp,categories(Yt),cmap,gmap,4);
save.pdf(sprintf('%s%sfs_%s_RVtest.ai',outDir,filesep,name),12,12);
plot.confusiongrad(YT,yc,cellstr(unique(Yt)),sprintf('df=%i',df),600,600,0.05,2,[0 0 0 1]);
save.pdf(sprintf('%s%sfs_%s_CMtest.ai',outDir,filesep,name),6,6);

[yc,yp] =ml.predict.glmnet(mdl,xxT,[],mdl.lambda(idx));
plot.radviz([yT,yc],yp,categories(Yt),cmap,cmap,4);
save.pdf(sprintf('%s%sfs_%s_RVmutest.ai',outDir,filesep,name),12,12);
plot.confusiongrad(yT,yc,cellstr(unique(Yt)),sprintf('df=%i',df),600,600,0.05,2,[0 0 0 1]);
save.pdf(sprintf('%s%sfs_%s_CMmutest.ai',outDir,filesep,name),6,6);

cmap = [151,84,204;129,198,121;243,184,99]./256;% frac/Num/Pos

tmap = hsv(360);
fmap = get.cmap([tmap(220,:);tmap(220,:).^.5;.85,.85,.85;cmap(1,:);cmap(1,:).^2]);
nmap = get.cmap([tmap(220,:);tmap(220,:).^.5;.85,.85,.85;cmap(2,:);cmap(2,:).^3]);
pmap = get.cmap([tmap(220,:);tmap(220,:).^.5;.85,.85,.85;cmap(3,:);cmap(3,:).^3]);

T = atlas.toTable(data.roi.roi.L,G);
plot.projectedAdj(mdl.fs,T,struct('or','a','cmap',cmap));
save.pdf(sprintf('%s%sfs_%s_mdl_Axial.ai',outDir,filesep,name),6,6);
plot.projectedAdj(mdl.fs,T,struct('or','rs','cmap',cmap));
save.pdf(sprintf('%s%sfs_%s_mdl_Saggital.ai',outDir,filesep,name),6,6);
plot.projectedAdj(mdl.fs,T,struct('or','c','cmap',cmap));
save.pdf(sprintf('%s%sfs_%s_mdl_Coronal.ai',outDir,filesep,name),6,6);

plot.multiCAdj(T,mdl.fs,cmap,10,2);
save.pdf(sprintf('%s%sfs_%s_Circ.ai',outDir,filesep,name),6,6);

c3nl.copy(c3nl.select('name','*.tex'),'/Users/eyalsoreq/Dropbox/PhD/CURRENT/papers/WM/wm-decoding-project/supplamentryMaterial')

%% 

fs = any(mdl.fs,2);
tbl = [data.task1.PPI(:,1:6),array2table(data.task1.PPI{:,7:end})];
W = grpstats(tbl,{'domain','load'},{'mean'},'DataVars',tbl.Properties.VariableNames(7:end));
cmap =[230,230,254;0,110,220;0,10,118]./256;
plot.multiGraph(TT,mdl.fs,categorical(W.domain).*categorical(W.load),W{:,4:end},{unique(W.domain),unique(W.load)},cmap,1);
save.pdf(sprintf('%s%sfs_%s.ai',outDir,filesep,'domainbyload_force'),22,6);
