clear all; close all; clc

%% path

addpath(genpath('C:\Users\u0098084\Documents\MATLAB\Metabolic Energy'));
addpath(genpath('C:\Users\u0098084\Documents\MATLAB\Metabolic Energy Tijs'));

%% input

input.subject = 'Leila_Arnalsteen';
input.initials = 'LA';
input.session = '20171120';
input.model = ['C:\Users\u0098084\Documents\Projecten\ATSECW\Processing\' input.subject '\' input.session '\Models\gait2392_simbody_ATSECW_' input.initials '_scaled.osim'];
input.speed = {'SLOW','FAST','COMF','6MWD'};
input.xlsx = fullfile('C:\Users\u0098084\Documents\Projecten\ATSECW\Data', input.subject, input.session, [input.subject, '.xlsx']);


%% compute metabolic energy

for i=1:length(input.speed)
    trial = [input.initials, '_', input.session, '_walking_', input.speed{i}, '_08'];
    do.phase(i) = load(['C:\Users\u0098084\Documents\Projecten\ATSECW\Processing\' input.subject '\' input.session '\DynamicOptimization\' trial '_DO_EN']);
    [energy.phase(i)] = computeMetabolicEnergy_2003_2010_2016(input,do.phase(i));
end

%% normaliseren tijd

for i=1:length(input.speed)
    energy.phase(i).cumsum_2003 = max(cumsum(energy.phase(i).model_2003));
    energy.phase(i).time_2003 = energy.phase(i).cumsum_2003/(do.phase(i).time(2)-do.phase(i).time(1));
    energy.phase(i).cumsum_2010 = max(cumsum(energy.phase(i).model_2010));
    energy.phase(i).time_2010 = energy.phase(i).cumsum_2010/(do.phase(i).time(2)-do.phase(i).time(1));
    energy.phase(i).cumsum_2016 = max(cumsum(energy.phase(i).model_2016));
    energy.phase(i).time_2016 = energy.phase(i).cumsum_2016/(do.phase(i).time(2)-do.phase(i).time(1));
end

%% compute experimental energy

[estat, energy.phase(1).exp, energy.phase(2).exp, energy.phase(3).exp, energy.phase(4).exp] = read_oxycon_output(input);

%% plot

figure()
for i = 1:length(input.speed)
    plot(energy.phase(i).exp, energy.phase(i).cumsum_2003, 'bo');hold on
    plot(energy.phase(i).exp, energy.phase(i).cumsum_2010, 'ro');hold on
    plot(energy.phase(i).exp, energy.phase(i).cumsum_2016, 'ko');hold on
    legend('2003', '2010', '2016')
end


%% save

%% path

rmpath(genpath('C:\Users\u0098084\Documents\MATLAB\Metabolic Energy'));
rmpath(genpath('C:\Users\u0098084\Documents\MATLAB\Metabolic Energy Tijs'));
