function [transitionProba_cell,time_cell1, nodeValue_cell, time_cell2] = transprob(n_nodes,time_series)
% Generate transition probability (cell) and 
% corresponding node values (cell) with corresponding time stage cells

n_scenario = size(time_series,1);
n_timestep = size(time_series,2);
transitionNode = zeros(n_nodes,n_nodes,n_timestep-1);
transitionProb = zeros(n_nodes,n_nodes,n_timestep-1);
transitionProba_cell = cell(n_timestep-1,1);

% time stage cell
time_cell1 = cell(1,n_timestep-1);
time_cell2 = cell(1,n_timestep);
for t = 1:(n_timestep-1)
    time_cell1{t} = num2str(t+1);
    time_cell2{t} = num2str(t);
end
time_cell2{n_timestep} = num2str(n_timestep);

% make ranges
value_ranges = zeros(n_nodes+1,n_timestep); % ranges of net loads
nodeValue = zeros(n_timestep,n_nodes); % node values
sort_ts = sort(time_series,1); % sorted net load
value_ranges(1,:) = sort_ts(1,:);
for i = 1:n_nodes
    for t = 1:n_timestep
        value_ranges(i+1,t) = sort_ts(round(n_scenario/n_nodes*i),t);
        if i == 1
           nodeValue(t,i) = mean(sort_ts(1:round(n_scenario/n_nodes*i),t)); % take the average within the range
        else
           nodeValue(t,i) = mean(sort_ts(round(n_scenario/n_nodes*(i-1)):round(n_scenario/n_nodes*i),t));
        end
    end
end

% Building transition probability
for i=1:n_scenario
    transisionCount = zeros(n_nodes,n_timestep); % count the pv_values that in the range
    for t=1:n_timestep
        for j=1:n_nodes
            if j ~= 1
                if value_ranges(j,t) < time_series(i,t) && time_series(i,t) <= value_ranges(j+1,t) 
                    transisionCount(j,t) = transisionCount(j,t)+1; 
                end
            else
                if value_ranges(j,t) <= time_series(i,t) && time_series(i,t) <= value_ranges(j+1,t) 
                    transisionCount(j,t) = transisionCount(j,t)+1; 
                end
            end
        end
    end
    for t = 2:n_timestep
        for l = 1:n_nodes
            for q = 1:n_nodes
                if transisionCount(l,t-1) == 1 && transisionCount(q,t) == 1
                    transitionNode(l,q,t-1) = transitionNode(l,q,t-1) + 1;
                end
            end
        end
    end
end

for t = 1:n_timestep-1
    for j = 1:n_nodes
        for k = 1:n_nodes
            if sum(transitionNode(j,:,t)) == 0
                transitionProb(j,k,t) = 0;
            else
                transitionProb(j,k,t) = transitionNode(j,k,t)/sum(transitionNode(j,:,t));
            end
        end
    end
end

% For stage 1: Choose the middile value and prob
nodeValue(1,1) = nodeValue(1,floor(n_nodes/2));
nodeValue(1,2:end) = nodeValue(1,1);

%%
nodeValue_cell = cell(n_timestep,1);
% transProba_cell = cell(n_timestep-1,1);
% new_transProba_cell{1} = 1;
nodeValue_cell{1} = nodeValue(1,1);
for t = 2:n_timestep
   [unique_val, idx] = unique(nodeValue(t,:));
   nodeValue_cell{t} = zeros(1,length(unique_val));
   [prev_unique_val,prev_idx] = unique(nodeValue(t-1,:));
   transitionProba_cell{t-1} = zeros(length(prev_idx),length(idx)); 
   for i = 1:length(unique_val)
       n = idx(i);
       same_idx = nodeValue(t,:)==nodeValue(t,n);
       nodeValue_cell{t}(i) = sum(nodeValue(t,same_idx));
       for j = 1:length(prev_unique_val)
           m = prev_idx(j);
           prev_same_idx = nodeValue(t-1,:)==nodeValue(t-1,m);
           transitionProba_cell{t-1}(j,i) = sum(sum(transitionProb(prev_same_idx,same_idx,t-1)));
       end
   end
end
%% delete elements with zero probabilty of landing 
for t = 2:n_timestep
    curr_n = size(transitionProba_cell{t-1},2) ;
    delete_set = [];
    for j = 1:curr_n
       if sum(transitionProba_cell{t-1}(:,j)) == 0.0
           delete_set = [delete_set , j];
       end
    end
    for i = 1:length(delete_set)
        j = delete_set(i);
       transitionProba_cell{t-1}(:,j) = [];
       transitionProba_cell{t}(j,:) = [];
       nodeValue_cell{t}(j) = []; 
    end
end
end