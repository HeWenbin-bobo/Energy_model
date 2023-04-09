%%%%%%%%Generate prediction model of ductile-to-brittle depth of cut of Silicon%%%%%%%%%%
%% Initialization
clc;
clear;
close all;

%% Start of time
tic;

%% Add related meta-heuristic algorithms
addpath(genpath('.\algorithms'));

%% Meta-heuristic algorithm initialization
algorithm_name = 'miga';
switch lower(algorithm_name)
    case 'miga'
        algorithm = @MIGA_main;
    case 'edo'
        algorithm = @EDO_main;
    otherwise
        algorithm = @MIGA_main;
end

%% Experiment result
experiment_result = readtable('.\input_data\experiment_results.xlsx','sheet',1);
input_features = readtable('.\input_data\experiment_results.xlsx','sheet',2);

%% Display setting
fontsize = 16;
set(0,'DefaultAxesFontSize',fontsize);
set(0,'DefaultTextFontSize',fontsize);

%% Parameter setting
Parameters.Nominal_cutting_speed = 1000;
Parameters.Nose_radius = 500;
Parameters.Edge_radius = 100;
Parameters.Rake_angle = -35; % 带正负号
Parameters.Flank_angle = 10;
Parameters.Rake_angle_rad = Parameters.Rake_angle*pi/180;
Parameters.nums = 1000; % 切深分段数量
Parameters.max_detph = 1E-6; % 单位m
Parameters.d_interval = Parameters.max_detph/Parameters.nums; % 切深分辨率
Parameters.plot_nums = fix(0.5E-6/Parameters.d_interval); % 绘图绘制多少个切深分辨率，fix是向零取整，总共绘制0.5E-6m即500nm切深
defaultParameters = Parameters;

%% Found truth value of Ft using meta-heuristic algorithm
DOC_truth = NaN(size(experiment_result,1),1);
DOC_pred = NaN(size(experiment_result,1),1);
fitness = NaN(size(experiment_result,1),1);
Ft_truth = NaN(size(experiment_result,1),1);
parfor i=1:size(experiment_result,1)
    Parameters = defaultParameters;
    % Parameter re-set
    Parameters.Rake_angle = experiment_result.rake_angle(i); % 带正负号
    Parameters.Rake_angle_rad = Parameters.Rake_angle*pi/180;
    Parameters.Nominal_cutting_speed = experiment_result.cutting_speed(i);
    Parameters.Laser_temperature = @calculate_net_equation_LT;
    Parameters.Laser_power = experiment_result.laser_power(4);
    Parameters.Crystal_direction = experiment_result.crystal_direction(i);
    Parameters.Temperature = Parameters.Laser_temperature(Parameters.Laser_power);

    Hardness = [0,12;400,11;600,3;750,1.5;900,0.75];
    Elastic_modulus = @(T) (129.09E9-0.01413E9*T*exp(-709/T))/1E9; % GPa
    Fracture_toughness = [25,0.8;900,0.75];
    Parameters.Hardness = linear_interpolation(Hardness,Parameters.Temperature);
    Parameters.Elastic_modulus = Elastic_modulus(Parameters.Temperature);
    Parameters.Fracture_toughness = linear_interpolation(Fracture_toughness,Parameters.Temperature);
    % Fracture_toughness与晶向应该有关
    Parameters.Surface_energy = 1.8;

    Parameters.v = 0.215; % 工件材料的泊松比，周行论文为v
    % Ft = 0.04; % 裂纹产生的临界载荷，周行论文为Ft
    Parameters.sigmay = 165E6; % 硅的屈服强度(Pa)

    % Calculate truth value of DOC (To be confirmed)
    temp = [experiment_result.DOC_1(i), experiment_result.DOC_2(i), experiment_result.DOC_3(i)];
    target = mean(temp(~isnan(temp(:,1)))); % 脆塑转变深度的目标值
    DOC_truth(i) = target; % Record truth value of DOC

    % Find minima of a function with meta-heuristic algorithm
    fobj = @(Ft) fun(Parameters, Ft, target); % Fitness function
    lb = 0; % Lower bound of Ft
    ub = 3; % Upper bound of Ft
    dim = 1; % Dimension of fobj input
    [ Ft_truth(i), fitness(i) ] = algorithm(lb,ub,dim,fobj); % Execute Algorithm

end

% Mark the end of find truth value of Ft
total_time = toc;

%% Genete a proxy model for predicting value of Ft using fitrauto
% rake_angle, cutting_speed, laser_power, crystal_direction为自变量
% Ft为因变量
variablename = {'rake_angle', 'cutting_speed', 'laser_power', 'crystal_direction', 'Ft'}; % 表头名字

% Generate training dataset
Ft = [experiment_result.order, Ft_truth];
index = Ft(:,1)~=-1;
Ft = Ft(index,2);
rake_angle = experiment_result.rake_angle(index);
cutting_speed = experiment_result.cutting_speed(index);
laser_power = experiment_result.laser_power(index);
crystal_direction = experiment_result.crystal_direction(index);

% 组成一个train table
Tbl_train = table(rake_angle, cutting_speed, laser_power, crystal_direction, Ft, 'VariableNames',variablename);

% Execute fitrauto to find a suitable proxy model automatically
% options = struct("Optimizer","asha","UseParallel",true,'Kfold',20); % "Kfold",5 if do not specify any cross-validation field
% options = struct("Optimizer","bayesopt","UseParallel",true,'Kfold',size(Tbl_train,1),"Repartition",true,"MaxObjectiveEvaluations",500); % "Kfold",5 if do not specify any cross-validation field
options = struct("Optimizer","bayesopt","UseParallel",true,'Kfold',size(Tbl_train,1),"Repartition",true); % "Kfold",5 if do not specify any cross-validation field
[Mdl,OptimizationResults] = fitrauto(Tbl_train,'Ft', ...
    'OptimizeHyperparameters','all','Learners','all', ...
    "HyperparameterOptimizationOptions",options, ...
    "CategoricalPredictors","crystal_direction");
% [Mdl,OptimizationResults] = fitrauto(Tbl_train,'Ft', ...
%     'OptimizeHyperparameters','all','Learners',["ensemble","kernel","linear","svm","net","tree"], ...
%     "HyperparameterOptimizationOptions",options, ...
%     "CategoricalPredictors","crystal_direction");
% [Mdl,OptimizationResults] = fitrauto(Tbl_train,'Ft', ...
%     'OptimizeHyperparameters','all','Learners',["ensemble","kernel","linear","svm","tree"], ...
%     "HyperparameterOptimizationOptions",options, ...
%     "CategoricalPredictors","crystal_direction");

% Generate test dataset
Ft = [experiment_result.order,Ft_truth];
index = Ft(:,1)==-1;
Ft = Ft(index,2);
rake_angle = experiment_result.rake_angle(index);
cutting_speed = experiment_result.cutting_speed(index);
laser_power = experiment_result.laser_power(index);
crystal_direction = experiment_result.crystal_direction(index);

% 组成一个test table
Tbl_test = table(rake_angle, cutting_speed, laser_power, crystal_direction, Ft, 'VariableNames',variablename);

% Calculate Ft prediction and error. Mdl is the proxy model.
Ft_pred = predict(Mdl,[experiment_result.rake_angle, experiment_result.cutting_speed, experiment_result.laser_power, experiment_result.crystal_direction]);
Ft_error = Ft_pred-Ft_truth;

%% Calculate DOC prediction according to Ft prediction
parfor i=1:size(experiment_result,1)
    Parameters = defaultParameters;
    % Parameter re-set
    Parameters.Rake_angle = experiment_result.rake_angle(i); % 带正负号
    Parameters.Rake_angle_rad = Parameters.Rake_angle*pi/180;
    Parameters.Nominal_cutting_speed = experiment_result.cutting_speed(i);
    Parameters.Laser_temperature = @calculate_net_equation_LT;
    Parameters.Laser_power = experiment_result.laser_power(4);
    Parameters.Crystal_direction = experiment_result.crystal_direction(i);
    Parameters.Temperature = Parameters.Laser_temperature(Parameters.Laser_power);

    Hardness = [0,12;400,11;600,3;750,1.5;900,0.75];
    Elastic_modulus = @(T) (129.09E9-0.01413E9*T*exp(-709/T))/1E9; % GPa
    Fracture_toughness = [25,0.8;900,0.75];
    Parameters.Hardness = linear_interpolation(Hardness,Parameters.Temperature);
    Parameters.Elastic_modulus = Elastic_modulus(Parameters.Temperature);
    Parameters.Fracture_toughness = linear_interpolation(Fracture_toughness,Parameters.Temperature);
    % Fracture_toughness与晶向应该有关
    Parameters.Surface_energy = 1.8;

    Parameters.v = 0.215; % 工件材料的泊松比，周行论文为v
    % Ft = 0.04; % 裂纹产生的临界载荷，周行论文为Ft
    Parameters.sigmay = 165E6; % 硅的屈服强度(Pa)

    % Calculate DOC prediction
    DOC_pred(i) = energy_model(Parameters, Ft_pred(i));

end

% Calculate DOC error
DOC_error = DOC_pred-DOC_truth;

%% Assemble all results
result = array2table([experiment_result.order, experiment_result.rake_angle, experiment_result.cutting_speed, experiment_result.laser_power, experiment_result.crystal_direction, ...
    DOC_truth, DOC_pred, DOC_error, Ft_truth, Ft_pred, Ft_error],...
    'VariableNames', {'order', 'rake_angle', 'cutting_speed', 'laser_power', 'crystal_direction', ...
    'DOC_truth', 'DOC_pred', 'DOC_error', 'Ft_truth', 'Ft_pred', 'Ft_error'});

save result;

%------------------------------------------%
function temperature = calculate_net_equation_LT(laser_power)
x=0.2*laser_power-1;
a1=2/(1+exp(-2*(78.3359*x-55.1104)))-1;
a2=2/(1+exp(-2*(-1.1104*x+1.1501)))-1;
a3=2/(1+exp(-2*(-77.6220*x-15.7035)))-1;
a4=2/(1+exp(-2*(1.5594*x+1.0230)))-1;
temperature=633.375*(a1*0.0128+a2*(-1.8417)+a3*(-0.005)+a4*0.1587+0.8985+1)+20;
end

%------------------------------------------%
function y = linear_interpolation(parameters, x)
last = 0;
next = 0;
for i=1:size(parameters,1)
    if x<parameters(i,1)
        next = i;
        break
    elseif x>parameters(i,1)
        last = i;
    else
        y = parameters(i,2);
        return
    end
end
if next==1
    y = parameters(1,2);
elseif last==size(parameters,1)
    y = parameters(end,2);
else
    x1 = parameters(last,1);
    x2 = parameters(last+1,1);
    y1 = parameters(last,2);
    y2 = parameters(last+1,2);
    a = (y2-y1)/(x2-x1);
    b = y1-a*x1;
    y = a*x+b;
end

end