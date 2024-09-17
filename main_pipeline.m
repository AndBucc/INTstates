%%% The topographic dynamic INT states analysis pipeline - by Andrea
%%% Buccellato. Remember that you can run separately sections starting with
%%% "%%" with the shortcut "CTRL+ENTER" (Win OS).

% !!! your dynamic INT matrices should have dimensions
% subjectsXchannelsXwindows (if obtained with sliding windows, otherwise
% it's timepoints instead of windows. See the slidingWindow() function included in this repo).

input = yourmatrix; % the input dynACW0 matrix. It should be a subjXchannelsXwindows cell/matrix
n_chans = size(input,2);

%% reshape input matrix according in preparation to k-means training - and delete zero columns

feed = zeros(n_chans,1);
for i = 1:size(input,1)
    temp = squeeze(input(i,:,:)); %subjXchanXwin
    feed = [feed,temp];
end


colzero = any(feed == 0, 1);
feed(:,colzero) = [];


% feed = (feed-mean(feed,2))./std(feed,0,2); % mean centering normalization
feed = (feed - min(feed)) ./ (max(feed) - min(feed)); % min-max normalization

%% k-means to find INT maps

% find best k (number of maps) as the elbow point of the k vs within-cluster
% variance graph
all_variances = cell(25,1);
all_Cs = cell(25,1);

for iter = 1:25

    for k = 2:20
        [idx{k},C{k},sumd{k},d{k}] = kmeans(feed',k,'Distance','cityblock','MaxIter',1000,'Replicates',1000);
    end

    all_variances{iter} = sumd;
    all_Cs{iter} = C;

end

%% obtain subject-wise histograms with back-fitting

% back-fitting

allgroup_maptimeseries = cell(size(input(1:end),1),1); %for all groups back-fitting

insp = {}; % not to delete - correlation matrices
% for backfitting from whole group training to single groups

for su = 1:size(input,1)
    one = squeeze(input(su,:,:));
    col_zero = any(one == 0, 1);
    one(:,col_zero) = []; % keep zero columns (windows) out
    if isempty(one)
        disp(['subject', num2str(su),'empty. Skipped'])
    else
        clear corre_int list_maps

        for i = 1:size(one,2)
            temp = one(:,i);
            for l = 1:size(C,1)
                %corre_int(l) = spatial_corr(C(l,:)', temp); %spatial correlation a là Murray. Same when data is mean detrended
                [corre_int(l),~]  = corr(C(l,:)', temp); %simple Pearson correlation

            end
            if all(isnan(corre_int))
                list_maps(i) = 0;
            else

                inspi(i,:) = corre_int;

                f = find(corre_int == max(corre_int));
                list_maps(i) = f;
            end

        end

        insp{su} = inspi;
        clear inspi

        %figure
        %histogram(list_maps)
        allgroup_maptimeseries_{su} = list_maps;
    end

end
%% characterize microstates
% ONLY AFTER MAPS ARE ESTABLISHED!!! OTHERWISE KEEP COMMENTED. Substitute
% "allgroup_maptimeseries" whenever needed at the beginning

% group-wise ApEn entropy
clear all_all
all_all = allgroup_maptimeseries(~cellfun('isempty',allgroup_maptimeseries)); % select only non-zero arrays in cell, just in case

%% Approximate Entropy (subj-wise) - parametrize as you like

m = 3;
clear ApEnsubj
for i = 1:size(all_all,2)
    temp = all_all{1,i};
    [App,phi] = ApEn(temp, 'm', 3, 'tau', 1);
    ApEnsubj(subj) = App(end);
end

disp(['avg ApEn:' num2str(mean(ApEnsubj))])

%% MEG only: Core-Periphery analysis


for map = 1:size(C,1)
    one = C(map,:);

    % periphery
    for i = 1:length(ito_periphery)
        indiana = find(ito == ito_periphery(i));
        periphery{i} = one(indiana);
    end

    %core
    for i = 1:length(ito_core)
        indiana = find(ito == ito_core(i));
        core{i} = one(indiana);
    end

    % compute C/P ratio
    newcore = [];
    for i = 1:8
        temp = core{1,i};
        newcore = [newcore,temp];
    end
    newperi = [];
    for i = 1:4
        temp = periphery{1,i};
        newperi = [newperi,temp];
    end
    ratios(map) = mean(newcore)/mean(newperi)
    clear newperi newcore
end