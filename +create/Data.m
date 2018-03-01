%% Get Behavioural data 

D = new.folder('Data/behaviour',pwd,{'results'});
system(sprintf( 'wget -O %s/result.mat %s',D.results,'https://www.dropbox.com/s/38z4owcboifeq1f/results.mat?dl=0'))

%% Move to copy over the conjunction files for domain general and stage general 
D = new.folder('Data/conjunction',pwd,{'results'});

clear
root = '/Users/eyalsoreq/Dropbox/PhD/CURRENT/papers/WM/wm-decoding-project/results/Data';
outDir = [root '/conjunction'];
% copy data from the cluster
c3nl.copy(c3nl.select('pth','/Users/eyalsoreq/mount/oldeyal/data/WM/paper/TASK1/H/h01','name','*Stage*_corrected.nii'),[outDir,'/task1']);
c3nl.copy(c3nl.select('pth','/Users/eyalsoreq/mount/oldeyal/data/WM/paper/TASK1/H/h01','name','*Stage*.csv'),[outDir,'/task1']);
c3nl.copy(c3nl.select('pth','/Users/eyalsoreq/mount/oldeyal/data/WM/paper/TASK1/H/h01','name','*Type*_corrected.nii'),[outDir,'/task1']);
c3nl.copy(c3nl.select('pth','/Users/eyalsoreq/mount/oldeyal/data/WM/paper/TASK1/H/h01','name','*Type*.csv'),[outDir,'/task1']);

c3nl.copy(c3nl.select('pth','/Users/eyalsoreq/mount/oldeyal/data/WM/paper/TASK2/H/h01','name','*Stage*_corrected.nii'),[outDir,'/task2']);
c3nl.copy(c3nl.select('pth','/Users/eyalsoreq/mount/oldeyal/data/WM/paper/TASK2/H/h01','name','*Stage*.csv'),[outDir,'/task2']);
c3nl.copy(c3nl.select('pth','/Users/eyalsoreq/mount/oldeyal/data/WM/paper/TASK2/H/h01','name','*Type*_corrected.nii'),[outDir,'/task2']);
c3nl.copy(c3nl.select('pth','/Users/eyalsoreq/mount/oldeyal/data/WM/paper/TASK2/H/h01','name','*Type*.csv'),[outDir,'/task2']);

%% Copy over mined data
outDir = [root '/MinedData'];
ls /Users/eyalsoreq/Dropbox/PhD/CURRENT/papers/WM/NN/Data/Updateddata/
c3nl.copy(c3nl.select('pth','/Users/eyalsoreq/Dropbox/PhD/CURRENT/papers/WM/NN/Data/Updateddata','name','[sSD]*.mat'),outDir);
% tidy the names etc...
fn = c3nl.select('pth',outDir,'name','*.mat');

for ii=1:numel(fn)
    load(fn{ii});
    data.task1.BA.Properties.VariableNames(2:6) = {'run','domain','load','state','rep'};
    data.task2.BA.Properties.VariableNames(2:7) = {'run','mainpulation','domain','load','state','rep'};
    data.task1.PPI.Properties.VariableNames(2:6) = {'run','domain','load','state','rep'};
    data.task2.PPI.Properties.VariableNames(2:7) = {'run','mainpulation','domain','load','state','rep'};
    save(fn{ii},'data','-V7.3');
end



%% decoding results 



%% Multivariate coding of WM factors for WB,DG and SG sets to be run on a cluster due to huge number of models run on the cluster
cd /group/hampshire_hub/oldeyal/data/WM/paper/ML
input = '/group/hampshire_hub/oldeyal/data/WM/paper/data1/';
sets = {'DG','SG','shen268'};
factors = {'domain','state','load'};
metrics = {'BA','PPI'};
for s = sets
    load ([input,s{1},'_gpp.mat']);
    for f = factors
        for m = metrics
            if c3nl.strDetect(m,'BA');idx = 2;else; idx = 1;end
            ix = ~(data.task1.(m{1}).state=='Cue');
            X1 = grpstats(data.task1.(m{1})(ix,:),{'subj',f{1},'run'},{'mean'},'DataVars',data.task1.(m{1}).Properties.VariableNames(6+idx:end));%
            ix = ~(data.task2.(m{1}).state=='Cue');
            X2 = grpstats(data.task2.(m{1})(ix,:),{'subj',f{1},'run'},{'mean'},'DataVars',data.task2.(m{1}).Properties.VariableNames((7+idx):end));%
            clear St
            St.name = sprintf('%s_%s_%s',s{1},f{1},m{1});
            St.X = X1{:,5:end};St.Y = categorical(cellstr(X1.(f{1})));
            St.stratify = X1.subj;
            St.XR = X2{:,5:end};St.YR = categorical(cellstr(X2.(f{1})));
            St.mc = 100;
            St.perm =10;
            run.ml.ecoc([],St);
        end
    end
end

%% Multivariate coding of WM domains at trial level with load as subgroups 

cd /group/hampshire_hub/oldeyal/data/WM/paper/ML
input = '/group/hampshire_hub/oldeyal/data/WM/paper/data1/';
sets = {'DG','SG','shen268'};
f = {'domain'};
metrics = {'BA','PPI'};
for s = sets
    load ([input,s{1},'_gpp.mat']);
    for m = metrics
        if c3nl.strDetect(m,'BA');idx = 2;else; idx = 1;end
        ix1 = ~(data.task1.(m{1}).state=='Cue');
        ix2 = ~(data.task2.(m{1}).state=='Cue');
        clear St
        St.name = sprintf('%s_%s_%s',s{1},'domainByLoad',m{1});
        X1 = grpstats(data.task1.(m{1})(ix1,:),{'subj','load','domain','run','rep'},{'mean'},'DataVars',data.task1.(m{1}).Properties.VariableNames(6+idx:end));%;
        St.X = X1{:,7:end};
        St.Y = categorical(cellstr(X1.domain));
        St.stratify = X1.subj;
        X2 = grpstats(data.task2.(m{1})(ix2,:),{'subj','load','domain','mainpulation','run','rep'},{'mean'},'DataVars',data.task2.(m{1}).Properties.VariableNames(7+idx:end));%;
        St.XR = X2{:,8:end};
        St.YR =  categorical(cellstr(X2.domain));
        St.mc = 100;
        St.perm =10;
        St.subgroup = struct('S',X1.load,'SR',X2.load);
        run.ml.ecoc([],St);
    end
end


%% Multivariate coding of WM domains over processing stages for all sets and metrics
sets = {'DG','SG','shen268'};
factors = {'domain','state','load'};
metrics = {'BA','PPI'};
for s = sets
    load ([input,s{1},'_gpp.mat']);
    for m = metrics
        if c3nl.strDetect(m,'BA');idx = 2;else; idx = 1;end
        ix = ~(data.task1.(m{1}).state=='Cue');
        X1 = grpstats(data.task1.(m{1})(ix,:),{'subj','domain','state','run'},{'mean'},'DataVars',data.task1.(m{1}).Properties.VariableNames(6+idx:end));%
        ix = ~(data.task2.(m{1}).state=='Cue');
        X2 = grpstats(data.task2.(m{1})(ix,:),{'subj','domain','state','run'},{'mean'},'DataVars',data.task2.(m{1}).Properties.VariableNames((7+idx):end));%
        clear St
        St.name = sprintf('%s_%s_%s',s{1},'DomainByState',m{1});
        St.X = X1{:,6:end};St.Y = categorical(cellstr(X1.domain));
        St.stratify = X1.subj;
        St.XR = X2{:,6:end};St.YR = categorical(cellstr(X2.domain));
        St.mc = 100;
        St.perm =10;   
        St.subgroup = struct('S',X1.state,'SR',X2.state);
        run.ml.ecoc([],St);
    end
end

%% 
%% Multivariate coding of mainpulation vs Maintain in study 2 
for s = sets
    load ([input,s{1},'_gpp.mat']);
    for m = metrics
        if c3nl.strDetect(m,'BA');idx = 2;else; idx = 1;end
        ix = data.task2.(m{1}).state=='Maintain';
        X = grpstats(data.task2.(m{1})(ix,:),{'subj','mainpulation','run'},{'mean'},'DataVars',data.task2.(m{1}).Properties.VariableNames((7+idx):end));%
        clear St
        St.name = sprintf('%s_%s_%s',s{1},'mainpulation',m{1});
        St.X = X{:,5:end};St.Y = categorical(cellstr(X.mainpulation));
        St.stratify = X.subj;
        St.mc = 1;
        St.perm =1;   
        run.ml.ecoc([],St);
    end
end


%% go over the decoding results - extract performance and create an avreged model

models = c3nl.select('pth','/group/hampshire_hub/oldeyal/data/WM/paper/ML/results','type','d','maxdepth',1,'mindepth',1);
T_factors = table();
T_domainByStages = table();
for m=models'
    fn = c3nl.select('pth',m{1},'name','*.mat');
    fp = strsplit(m{1},{'_','/'});
    for ii = 1:numel(fn)
        tmp = load(fn{ii});
        fp = strsplit(fn{ii},{'_','/','.'});
        if c3nl.strDetect(m,'By')
            T_domainByStages =[T_domainByStages;[cell2table(repmat(fp(end-5:end-2),size(tmp.output.SG,1),1),'VariableNames',{'perm','set','factor','metric'}),tmp.output.SG]];  
        else
            T_factors = [T_factors;[cell2table(repmat(fp(end-5:end-2),size(tmp.output.F1,1),1),'VariableNames',{'perm','set','factor','metric'}),...
                array2table(tmp.output.F1,'VariableNames',{'Null','oob','test','replication'})]];   
        end
        fprintf('*');
    end
end


T_factors = table();
models = c3nl.select('pth','/group/hampshire_hub/oldeyal/data/WM/paper/ML/results','type','d','maxdepth',1,'mindepth',1);
models = models(c3nl.strDetect(models,'ByLoad'));
for m=models'
    fn = c3nl.select('pth',m{1},'name','*.mat');
    fp = strsplit(m{1},{'_','/'});
    for ii = 1:numel(fn)
        tmp = load(fn{ii});
        fp = strsplit(fn{ii},{'_','/','.'});
        T_factors = [T_factors;[cell2table(repmat(fp(end-5:end-2),size(tmp.output.F1,1),1),'VariableNames',{'perm','set','factor','metric'}),...
                array2table(tmp.output.F1,'VariableNames',{'Null','oob','test','replication'})]]; 
        
    fprintf('*');
    end
end



M_perf = table();
models = c3nl.select('pth','/group/hampshire_hub/oldeyal/data/WM/paper/ML/results','type','d','maxdepth',1,'mindepth',1);
models = models(c3nl.strDetect(models,'ByLoad'));
for m=models'
    fn = c3nl.select('pth',m{1},'name','*.mat');
    fp = strsplit(m{1},{'_','/'});
    for ii = 1:numel(fn)
        t = load(fn{ii});
        fp = strsplit(fn{ii},{'_','/','.'});
        ix =  c3nl.strDetect(t.output.SG.model,'global');
        null =  t.output.F1(:,1);
        tmp = t.output.SG(ix,:);
        ixl = c3nl.strDetect(tmp.subgroup,'Low');
        ixm = c3nl.strDetect(tmp.subgroup,'Med');
        ixh = c3nl.strDetect(tmp.subgroup,'High');
        M_perf = [M_perf;[cell2table(repmat(fp(end-5:end-2),size(t.output.F1,1),1),'VariableNames',{'perm','set','factor','metric'}),...
                array2table([null,tmp.F1(ixl),tmp.F1(ixm),tmp.F1(ixh),tmp.F1_R(ixl),tmp.F1_R(ixh)],'VariableNames',{'Null','Low_cv','Med_cv','High_cv','Low_R','High_R'})]];   
        fprintf('*');

    end
end

save('loadPref2.mat','M_perf','T_factors')


