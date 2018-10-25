
clear all; close all; clc

%% path

addpath(genpath('C:\Users\u0098084\Documents\MATLAB\Metabolic Energy'));

%% load graphs

graph(1).name = 'graph2A';
graph(2).name = 'graph2B';
graph(3).name = 'graph2C';
graph(4).name = 'graph2D';

for i = 1:length(graph)
    graph(i).pattern = xlsread(['C:\Users\u0098084\Documents\MATLAB\Metabolic Energy\MetabolicEnergy\Holt2014\' graph(i).name '_pattern.csv']);
    graph(i).pattern(1,:) = 0;
    graph(i).pattern(end,:) = [0.4,0];
    graph(i).time = (0.001:(0.4-0)/400:0.4)';
    graph(i).strain = interp1(graph(i).pattern(:,1),graph(i).pattern(:,2),graph(i).time);
end

figure()
for i = 1:length(graph)
    subplot(length(graph),1,i)
    plot(graph(i).time,graph(i).strain);
    xlim([0 0.4]);
    ylim([-0.1 0.1])
end

for i = 1:length(graph)
    graph(i).exc = xlsread(['C:\Users\u0098084\Documents\MATLAB\Metabolic Energy\MetabolicEnergy\Holt2014\' graph(i).name '_exc.csv']);
    graph(i).exc_1 = zeros(size(graph(i).strain));
    graph(i).exc_1(find(graph(i).time>graph(i).exc(1),1,'first'):find(graph(i).time<graph(i).exc(2),1,'last'),1) = 1;
    graph(i).exc_1(find(graph(i).time>graph(i).exc(3),1,'first'):find(graph(i).time<graph(i).exc(4),1,'last'),1) = 1;
end

for i = 1:length(graph)
    subplot(length(graph),1,i)
    plot(graph(i).time,graph(i).strain, 'k'); hold on
    plot(graph(i).time(graph(i).exc_1==1,:),graph(i).strain(graph(i).exc_1==1,:),'ok', 'MarkerFaceColor', 'k');
    xlim([0 0.4]);
    ylim([-0.1 0.1])
end

%% muscle parameters

params(6,1) = 0.0152; % musclemass(g) - defined on page 4366 'results' 
params(7,1) = 1.056; % muscledensity(g/cm3) - defined in Kargo 2002 functional morphology
params(8,1) = 42.6; % maximal stress(N/cm2) - defined on page 4366 'results'
params(9,1) = 0; % tendons stiffness (0 to ignore the influence of the tendon)

params(2,1) = 0.0135; % lmopt(m) - defined on page 4366 'results' % what to do with fiber length = 0.93 lmopt?
params(1,1) = params(8,1)*(params(6,1)/params(7,1))/(params(2,1)*100); % maximal isometric force
params(3,1) = 0; % tendonslacklength - currently not used
params(4,1) = 0; % pennation angle(°) - currently not used
params(5,1) = 5.8*params(2,1); % maximal contraction velocity (m/s)
params(10,1)= 0.5; % pct slow twitch - difficult to find exact numbers - "typically mixed muscle fiber types in frog iliofibularis" - Muscle-nerve interaction book - Gerta Vrbova 

%% calculate activations

td = 0.06;
ta = 0.015;
a0 = 0;

for i = 1:length(graph)
    [graph(i).t,graph(i).a] = calculate_activations_excitations(graph(i).time, graph(i).exc_1,ta,td,a0);
    graph(i).ai = interp1(graph(i).t,graph(i).a,graph(i).time);
end

figure()
for i = 1:length(graph)
    subplot(length(graph),1,i)
    plot(graph(i).time,graph(i).ai); hold on
end

%% calculate muscle forces & muscle stress

for i = 1:length(graph)
    [graph(i).lmtilde, graph(i).vmtilde] = get_lmtilde_vmtilde(graph(i).strain,graph(i).time);
    graph(i).lmt = graph(i).lmtilde.*params(2,1);
    
%     figure()
%     subplot(211)
%     plot(graph(i).lmtilde)
%         
%     subplot(212)
%     plot(graph(i).vmtilde)

    % Parameters of active muscle force-velocity characteristic
    load ActiveFVParameters.mat
    Fvparam(1) = 1.475*ActiveFVParameters(1);
    Fvparam(2) = 0.25*ActiveFVParameters(2);
    Fvparam(3) = ActiveFVParameters(3) + 0.75;
    Fvparam(4) = ActiveFVParameters(4) - 0.027;

    % Parameters of active muscle force-length characteristic
    load Faparam.mat                            

    % Parameters of passive muscle force-length characteristic
    e0 = 0.6; kpe = 4; t50 = exp(kpe * (0.2 - 0.10e1) / e0);
    pp1 = (t50 - 0.10e1); t7 = exp(kpe); pp2 = (t7 - 0.10e1);
    Fpparam = [pp1;pp2];
    [graph(i).fm,graph(i).ft,graph(i).err,graph(i).fce, graph(i).fpe] = ForceEquilibrium_lMtildeState_holt2014(graph(i).ai,graph(i).lmtilde,graph(i).vmtilde,graph(i).lmt,params,Fvparam,Fpparam,Faparam);
    graph(i).sigma = graph(i).fm./(params(1,1)/params(8,1));
end

figure()
for i = 1:length(graph)
    subplot(length(graph),1,i)
    plot(graph(i).time,graph(i).sigma); hold on
end

%% calculate heat

for i = 1:length(graph)
    graph(i).fmltilde = get_fmltilde_holt(graph(i).lmtilde,Faparam);
    b = 1000; % smoothness models
    [graph(i).heat_2003, graph(i).heat_2003_am, graph(i).heat_2003_sl] = calculateheat_2003(graph(i).exc_1,graph(i).ai,graph(i).lmtilde,graph(i).vmtilde,params(10,1),params(5,1),graph(i).fmltilde,b); % based on Umberger 2003, energy due to scaling factor for activation?
    graph(i).heat_2003_cum = cumsum(graph(i).heat_2003).*0.001; % convert heat produced from W/kg to J/kg (*0.001)
    [graph(i).heat_2010, graph(i).heat_2010_am, graph(i).heat_2010_sl] = calculateheat_2010(graph(i).exc_1,graph(i).ai,graph(i).lmtilde,graph(i).vmtilde,params(10,1),params(5,1),graph(i).fmltilde,b); % based on Umberger 2003, energy due to scaling factor for activation?
    graph(i).heat_2010_cum = cumsum(graph(i).heat_2010).*0.001; % convert heat produced from W/kg to J/kg (*0.001)
    [graph(i).heat_2016, graph(i).heat_2016_am, graph(i).heat_2016_sl] = calculateheat_2016(graph(i).exc_1,graph(i).ai,graph(i).lmtilde,graph(i).vmtilde,params(10,1),params(5,1),graph(i).fmltilde,b); % based on Umberger 2003, energy due to scaling factor for activation?
    graph(i).heat_2016_cum = cumsum(graph(i).heat_2016).*0.001; % convert heat produced from W/kg to J/kg (*0.001)
    [graph(i).heat_2004, graph(i).heat_2004_am, graph(i).heat_2004_sl] = calculateheat_2004(graph(i).exc_1,graph(i).lmtilde,graph(i).vmtilde,params(10,1),graph(i).fm);
    graph(i).heat_2004_cum = cumsum(graph(i).heat_2004).*0.001;
end

figure()
for i = 1:length(graph)
    subplot(length(graph),1,i)
    plot(graph(i).time,graph(i).heat_2003); hold on
    plot(graph(i).time,graph(i).heat_2004); hold on
    plot(graph(i).time,graph(i).heat_2010); hold on
    plot(graph(i).time,graph(i).heat_2016);
    legend('2003', '2004', '2010', '2016')
    title(['heat production ',graph(i).name])
end

figure()
for i = 1:length(graph)
    subplot(length(graph),1,i)
    plot(graph(i).time,graph(i).heat_2003_cum); hold on
    plot(graph(i).time,graph(i).heat_2004_cum); hold on
    plot(graph(i).time,graph(i).heat_2010_cum); hold on
    plot(graph(i).time,graph(i).heat_2016_cum);
    legend('2003', '2004', '2010', '2016')
    title(['cumulative heat production ',graph(i).name])
end

figure()
for i = 1:length(graph)
    subplot(length(graph),1,i)
    plot(graph(i).time,graph(i).heat_2003_am); hold on
    plot(graph(i).time,graph(i).heat_2004_am); hold on
    plot(graph(i).time,graph(i).heat_2010_am); hold on
    plot(graph(i).time,graph(i).heat_2016_am);
    legend('2003', '2004', '2010', '2016')
    title(['am heat production ',graph(i).name])
end

figure()
for i = 1:length(graph)
    subplot(length(graph),1,i)
    plot(graph(i).time,graph(i).heat_2003_sl); hold on
    plot(graph(i).time,graph(i).heat_2004_sl); hold on
    plot(graph(i).time,graph(i).heat_2010_sl); hold on
    plot(graph(i).time,graph(i).heat_2016_sl);
    legend('2003', '2004', '2010', '2016')
    title(['sl heat production ',graph(i).name])
end

%% calculate net work

for i = 1:length(graph)
    graph(i).work = (graph(i).fm.*graph(i).vmtilde.*params(2,1))./params(6,1);
end

%% plot

figure()
for i = 1:length(graph)
    if i == 1; start = 1;
    elseif i == 2; start = 2;
    elseif i == 3; start = 7;
    elseif i == 4; start = 8;
    end
    subplot(length(graph)/2*3,length(graph)/2,start)
    plot(graph(i).time,graph(i).strain, 'k'); hold on
    plot(graph(i).time(graph(i).exc_1==1,:),graph(i).strain(graph(i).exc_1==1,:),'ok', 'MarkerFaceColor', 'k');
    xlim([0 0.4]);
    subplot(length(graph)/2*3,length(graph)/2,start+2)
    plot(graph(i).time,graph(i).sigma);
    xlim([0 0.4]);
    subplot(length(graph)/2*3,length(graph)/2,start+4)
    plot(graph(i).time,graph(i).heat_2003_cum); hold on
    plot(graph(i).time,graph(i).heat_2004_cum);
    plot(graph(i).time,graph(i).heat_2010_cum);
    plot(graph(i).time,graph(i).heat_2016_cum);
    legend('2003', '2004', '2010', '2016')
    xlim([0 0.4]);
end

% figure()
% for i = 1:length(graph)
%     graph(i).barplot(:,(i-1)*3+1) = sum(graph(i).work);
%     graph(i).barplot(:,(i-1)*3+2) = sum(graph(i).heat_2003);
%     graph(i).barplot(:,(i-1)*3+3) = graph(i).barplot(:,(i-1)*3+1)+graph(i).barplot(:,(i-1)*3+2);
%     bar(graph(i).barplot);
% end

%% path

rmpath(genpath('C:\Users\u0098084\Documents\MATLAB\Metabolic Energy'));
