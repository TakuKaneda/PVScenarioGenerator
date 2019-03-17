%% Create lattice with simulation data 
% use https://github.com/kyamagu/matlab-json/ for write .json file

n_nodes = 10;

% data of simulation
n_scenario = 5000;
N = 1000;
datafilename = 'PGE-SASH-4101_winter';
filename = ['../pv_simulation/',datafilename,'_',num2str(n_scenario),'sim_',num2str(N),'meshpdf.csv'];
pv_scenario = load(filename);

% normalize data
pv_capacity = 2.9;
pv_scenario_norm = pv_scenario/pv_capacity;

% build the lattice
[transitionProba_cell,time_cell1, nodeValue_cell, time_cell2] = transprob(n_nodes,pv_scenario_norm);
% save as .json files
TP_jsondata = jsonencode(containers.Map(time_cell1,transitionProba_cell));
Node_jsondata = jsonencode(containers.Map(time_cell2,nodeValue_cell));
json.write(TP_jsondata,['TP_',datafilename,'.json']);
json.write(Node_jsondata,['Node_',datafilename,'.json']);
%% node value in lattice
% figure()
% hold on
% timesteps = linspace(0,23.75,size(pv_scenario_norm,2));
% for t = 1:size(pv_scenario_norm,2)
%     num = length(nodeValue_cell{t});
%     scatter(timesteps(t).*ones(num,1),nodeValue_cell{t},'b','X');
%     line([timesteps(t) timesteps(t)], [nodeValue_cell{t}(1) nodeValue_cell{t}(end)]...
%         ,'Color','black','LineStyle','-');
% end
% xlabel('time (h)')
% ylabel('normalized power (kW)')
%% transition probability in lattice
% figure()
% for t = 2:size(pv_scenario_norm,2)
%     prev_num = size(transitionProba_cell{t-1},1);
%     num = size(transitionProba_cell{t-1},2);
%     for i = 1:prev_num
%         for j = 1:num
%             if transitionProba_cell{t-1}(i,j) > 0
%                line([timesteps(t-1) timesteps(t)],...
%                    [nodeValue_cell{t-1}(i) nodeValue_cell{t}(j)],...
%                    'LineWidth',transitionProba_cell{t-1}(i,j)*3)
%             end
%         end
%     end
% end
% xlabel('time (h)')
% ylabel('normalized power (kW)')