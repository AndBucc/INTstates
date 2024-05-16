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

    for i = 2:30
        [t_idx{i},t_C{i},sumd] = kmeans(feed', i, 'Distance', 'cityblock','MaxIter',2000); %here you can parametrize as preferred
        sums(i) = sum(sumd);
    end

    [res_x, idx_of_result] = knee_pt(sums(2:end), 2:30, true); %find knee point
    wethebestmusic(iter) = res_x; % store best k for this iteration

    all_variances{iter} = sums;
    all_Cs{iter} = t_C;

end
%% select k-means results with the best k


idx = t_idx{res_x};
C = t_C{res_x};

clear t_C t_idx

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
%% characterize microstates - entropy, distributions, etc..
% ONLY AFTER MAPS ARE ESTABLISHED!!! OTHERWISE KEEP COMMENTED. Substitute
% "allgroup_maptimeseries" whenever needed at the beginning

% group-wise average empirical and permutation entropy
clear all_all
all_all = allgroup_maptimeseries(~cellfun('isempty',allgroup_maptimeseries)); % select only non-zero arrays in cell, just in case

%% Permutation Entropy (subj-wise) - parametrize as you like

m = 5;
clear PermEnsubj
for i = 1:size(all_all_rest,2)
    temp = all_all_rest{1,i};
    [~, Pnorm, ~] = PermEn(temp, 'm', m, 'tau', 1, 'Norm', true);
    PermEnsubj(i) = Pnorm(end);
end

disp(['avg PermEn:' num2str(mean(PermEnsubj))])

%% Distribution Entropy (just for curiosity, not used in the study)
group_dist = [];
for i = 1:size(all_all,2)
    temp = all_all{1,i};
    group_dist = [group_dist,temp];
end
clear i

max_en = log(size(C,1));
[group_en, ~] =  DistEn(group_dist, 'm', 2, 'tau', 1, 'Bins', 6, 'Logx', 2, 'Norm',true);
disp(['Group Distribution Entropy:', num2str(group_en), newline, 'Maximal Entropy for n states:', num2str(max_en)])


%% average state duration

for i = 1:size(subj,1)

    a = subj(i,:);

    if isempty(a)
    end

    for j = 1:size(C,1)
        percentages(i,j) = (length(find(a == j))/size(input,3)) * 100;
    end
end

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

%% generate shuffled matrices for bootstrap statistics
all_shuffled = {};
for maps = 1:size(all_all,2)
    a = all_all{1,maps};
    n_rep = 1000;
    shuffled = zeros(n_rep, length(a));
    for rep = 1:n_rep
        indx = randperm(length(a));
        shuffled(rep,:) = a(indx);
    end
    all_shuffled{maps} = shuffled;
  
end