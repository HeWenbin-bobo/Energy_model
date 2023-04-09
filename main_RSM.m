%%%%%%%%The simulation of cutting%%%%%%%%%%
%% % https://ww2.mathworks.cn/help/stats/linearmodel.html
%%%%%Cutting parameters%%%%%
clc;
clear;
close all;

tic;

algorithm_name = 'miga';
switch lower(algorithm_name)
    case 'miga'
        algorithm = @MIGA_main;
    case 'edo'
        algorithm = @EDO_main;
    otherwise
        algorithm = @MIGA_main;
end

experiment_result = readtable('.\input_data\experiment_results.xlsx','sheet',1);
input_features = readtable('.\input_data\experiment_results.xlsx','sheet',2);

fontsize = 16;
set(0,'DefaultAxesFontSize',fontsize);
set(0,'DefaultTextFontSize',fontsize);

Parameters.Nominal_cutting_speed = 1000;
Parameters.Nose_radius = 500;
Parameters.Edge_radius = 100;
Parameters.Rake_angle = -35; % 带正负号
Parameters.Flank_angle = 10;
Parameters.Rake_angle_rad = Parameters.Rake_angle*pi/180;
Parameters.nums = 10000; % 切深分段数量
Parameters.max_detph = 1E-6; % 单位m
Parameters.d_interval = Parameters.max_detph/Parameters.nums; % 切深分辨率
Parameters.plot_nums = fix(0.5E-6/Parameters.d_interval); % 绘图绘制多少个切深分辨率，fix是向零取整

%%
bestind = NaN(size(experiment_result,1),7);
bestfit = NaN(size(experiment_result,1),2);
for i=1:size(experiment_result,1)
    %%
    Parameters.Rake_angle = experiment_result{i,2}; % 带正负号
    Parameters.Rake_angle_rad = Parameters.Rake_angle*pi/180;
    Parameters.Nominal_cutting_speed = experiment_result{i,3};
    Parameters.Laser_temperature = @calculate_net_equation_LT;
    Parameters.Laser_power = experiment_result{i,4};
    Parameters.Crystal_direction = experiment_result{i,5};
    Parameters.Temperature = Parameters.Laser_temperature(Parameters.Laser_power);

    %%
    Hardness = [0,12;400,11;600,3;750,1.5;900,0.75];
    Elastic_modulus = @(T) (129.09E9-0.01413E9*T*exp(-709/T))/1E9; % GPa
    Fracture_toughness = [25,0.8;900,0.75];
    Parameters.Hardness = linear_interpolation(Hardness,Parameters.Temperature);
    Parameters.Elastic_modulus = Elastic_modulus(Parameters.Temperature);
    Parameters.Fracture_toughness = linear_interpolation(Fracture_toughness,Parameters.Temperature);
    % Fracture_toughness与晶向应该有关
    Parameters.Surface_energy = 1.8;

    %%
    Parameters.v = 0.215; % 工件材料的泊松比，周行论文为v
    % Ft = 0.04; % 裂纹产生的临界载荷，周行论文为Ft
    Parameters.sigmay = 165E6; % 硅的屈服强度(Pa)

    %%
    temp = experiment_result{i,end-2:end};
    target = mean(temp(~isnan(temp))); % 脆塑转变深度的目标值

    %%
    bestind(i,1:5) = experiment_result{i,1:5};
    bestind(i,6) = target;
    bestfit(i,1) = experiment_result{i,1};

    %% Follow the example to implement MIGA， https://github.com/OpenLlop/MOT
    % Find minima of a function with Islands Model + Genetic Algorithm (GA)
    fobj = @(Ft) fun(Parameters, Ft, target); % Fitness function

    %
    lb = 0;
    ub = 3;
    dim = 1;

    % Execute Algorithm
    [ bestind(i,7), bestfit(i,2) ] = algorithm(lb,ub,dim,fobj);

end

total_time = toc;

%%
variablename = experiment_result.Properties.VariableNames(2:5); % 变量名提取
variablename = cellstr(variablename);
variablename(end+1) = {'Ft'};

% %% RSM拟合
% X = data(:,1:end-1);
% Y = data(:,end);
% rstool(X,Y,'quadratic'); % quadratic表明拟合公式为完全二次，数据只有16个，不够

%% Stepwise Regression+backward elimination 逐步回归+向后消除
% rake_angle, cutting_speed, laser_power, crystal_direction为自变量
% Ft为因变量
Ft = bestind(:,[1,7]);
index = Ft(:,1)~=-1;
Ft = Ft(index,2);
rake_angle = experiment_result.rake_angle(index);
cutting_speed = experiment_result.cutting_speed(index);
laser_power = experiment_result.laser_power(index);
crystal_direction = experiment_result.crystal_direction(index);

% 组成一个table
tbl = table(rake_angle, cutting_speed, laser_power, crystal_direction, Ft, 'VariableNames',variablename);
% 从二阶表达式出发进行向后消除，上限为二阶表达式（即得到的结果不会有三阶以上的项），因变量定义为Ft
% mdl1 = stepwiselm(tbl,'quadratic','PredictorVars',{'rake_angle','cutting_speed','laser_power','crystal_direction'},'CategoricalVar',{'crystal_direction'},'ResponseVar','Ft','Verbose',0);
mdl1 = stepwiselm(tbl,'quadratic','ResponseVar','Ft','Verbose',0);

disp(' ');

% 从三阶表达式出发进行向后消除，不设定上限模型，Verbose为2说明将优化过程给print出来
% mdl2 = stepwiselm(tbl,'poly3333','ResponseVar','Ft','Upper','poly3333','Verbose',2);
% mdl2 = stepwiselm(tbl,'poly3333','CategoricalVar',{'crystal_direction'},'ResponseVar','Ft','Verbose',0);
mdl2 = stepwiselm(tbl,'poly3333','ResponseVar','Ft','Verbose',0);


bestind(:,end+1) = predict(mdl2,bestind(:,2:5));
bestind(:,end+1) = bestind(:,end)-bestind(:,end-1);
bestind_size = size(bestind,2);

for i=1:size(experiment_result,1)
    %%
    Parameters.Rake_angle = experiment_result{i,2}; % 带正负号
    Parameters.Rake_angle_rad = Parameters.Rake_angle*pi/180;
    Parameters.Nominal_cutting_speed = experiment_result{i,3};
    Parameters.Laser_temperature = @calculate_net_equation_LT;
    Parameters.Laser_power = experiment_result{i,4};
    Parameters.Crystal_direction = experiment_result{i,5};
    Parameters.Temperature = Parameters.Laser_temperature(Parameters.Laser_power);

    %%
    Hardness = [0,12;400,11;600,3;750,1.5;900,0.75];
    Elastic_modulus = @(T) (129.09E9-0.01413E9*T*exp(-709/T))/1E9; % GPa
    Fracture_toughness = [25,0.8;900,0.75];
    Parameters.Hardness = linear_interpolation(Hardness,Parameters.Temperature);
    Parameters.Elastic_modulus = Elastic_modulus(Parameters.Temperature);
    Parameters.Fracture_toughness = linear_interpolation(Fracture_toughness,Parameters.Temperature);
    % Fracture_toughness与晶向应该有关
    Parameters.Surface_energy = 1.8;

    %%
    Parameters.v = 0.215; % 工件材料的泊松比，周行论文为v
    % Ft = 0.04; % 裂纹产生的临界载荷，周行论文为Ft
    Parameters.sigmay = 165E6; % 硅的屈服强度(Pa)

    bestind(i,bestind_size+1) = energy_model(Parameters, bestind(i,bestind_size-1));

end

bestind(:,end+1) = bestind(:,end)-bestind(:,end-4);

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