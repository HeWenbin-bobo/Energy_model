function critical_doc = energy_model(Parameters, Ft)

Rake_angle_rad = Parameters.Rake_angle_rad;

nums = Parameters.nums; % 切深分段数量
max_detph = Parameters.max_detph; % 单位m
d_interval = Parameters.d_interval; % 切深分辨率

rou = 2.1e6; % 硅的密度(kg/m^3)
Mphase = 28; % 硅的摩尔质量(g/mol)
Qc = 7e3; % 非晶硅的晶化能(J)
k2 = 10; % 度量常数（可调）
% d = linspace(0,max_detph,nums); %切深(m)
d = 0:d_interval:max_detph; %切深(m)
r = Parameters.Edge_radius.*10^-9; %刃口半径(m)
R = Parameters.Nose_radius.*10^-6; %圆弧半径(m)
vc = Parameters.Nominal_cutting_speed.*10^-3/60; %切削速度(m/s)
alpha = Parameters.Flank_angle*pi/180; %刀具后角(rad)
gamma = abs(Parameters.Rake_angle)*pi/180; %刀具前角(rad)
rc = 0.3; %Cutting ratio（未变形切屑厚度与切屑厚度的比值）
miu = 0.6; %刀具刃口与工件间摩擦系数（常规切削时）
beta = atan(miu); % 摩擦角（常规切削时）(rad)
H = Parameters.Hardness.*10^9; %硬度(Pa)
E = Parameters.Elastic_modulus.*10^9; %弹性模量(Pa)
Kc = Parameters.Fracture_toughness.*10^6; %断裂韧性(Pa)
Er = 0.6; %弹性变形率
k1 = 1.25; %后刀面与工件间界面应力的最佳比例常数
ys = Parameters.Surface_energy; %比表面能(J/m^2)
X = 0.064; %几何常数
C1 = 0.1458; %无量纲系数，与采用的何种材料和压头无关，用以计算裂纹相关指标
C2 = 0.1450; %无量纲系数，与采用的何种材料和压头无关，用以计算裂纹相关指标，论文是C1和C2用同一个值（周行这样写，估计是为了尝试不同值所带来的影响）
dp = r/2; %临界耕犁深度(m)
gammae = asin(1-d./r).*(d<dp)+ gamma.*(d>=dp);
b = (pi/2+gammae-alpha)/2; %刀具相对两侧边之间的半角
phie = (atan(rc.*cos(gammae)./(1-rc.*sin(gammae)))).*(d<dp) + (pi/4-beta/2+Rake_angle_rad/2).*(d>=dp); %等效剪切角(3-3)
s = (Er.*dp).*(d>dp)+(Er.*d).*(d<=dp); %弹性恢复量(3-7)
W = 2.*sqrt(dp.*(2*R-dp)).*(d>dp)+2.*sqrt(d.*(2*R-d)).*(d<=dp); %切削宽度(3-8)
Af = W.*(2.*s./3./sin(alpha)); %刀具后刀面的接触面积(3-6)
deltaf = k1.*H.*sqrt(H./E); % 挤压应力的经验关系，周行论文为deltaf(3-9)（要看这个挤压应力的方向是垂直于切削方向还是垂直于后刀面了，因为有些SCI论文里是按垂直于切削方向算的，方向不同，分力计算就不同）
miuf = (miu).*(d>=dp)+((cos(phie)-sin(phie)./sqrt(3))./(cos(phie)./sqrt(3)+sin(phie))).*(d<dp); %刀具工件之间的摩擦系数。因为切深小于tp也就是最小切削厚度的时候，是没有切屑的，所以利用剪切角来计算(3-10)。似乎应该是前后刀面共用一个摩擦系数了，常规切削的摩擦系数可以用实验来测
Ff = deltaf.*Af; % 后刀面形成法向正挤压力
theta = atan(sqrt(d.*(2.*R-d))./(R-d)); % 在任一切削深度处刀尖圆弧段在切削深度平面上对应的扇形角，周行论文为xita(3-13)
theta0 = atan(sqrt(dp.*(2.*R-dp))./(R-dp)); % 在临界耕犁深度处刀尖圆弧段在切削深度平面上对应的扇形角，周行论文为xita0(3-14)
Ac = R.^2.*theta-(R-d).*sqrt(d.*(2.*R-d)); % 工件材料瞬时去除量在切削方向投影面积(犁耕+剪切，也就是塑性)(3-11)
As = 0.*(d<dp)+(R.^2.*theta-(R-d).*sqrt(d.*(2.*R-d))-R.^2.*theta0+(R-dp).*sqrt(dp.*(2.*R-dp))).*(d>=dp); % 剪切主导区域未切削工件材料在切削方向投影面积(3-12)。因为切深较小的时候，只有犁耕，所以前面有些值是负数，将它们设置成0

%% 耕犁和塑性阶段
% 公式第一项
% 耕犁阶段
int_part = @(t) cot(atan(rc.*cos(asin(1-t./r))./(1-rc.*sin(asin(1-t./r))))).*2.*sqrt(t.*(2.*R-t)); % 积分项
Fcplunge = zeros(1, numel(d));
Fcplungec = H./3./sqrt(3).*integral(int_part, 0, dp); % 临界值
for i = 1:numel(d) % if cutting depth <= critical doc, then calculate by integral, otherwise equals to 临界值
    if d(i)<dp
        Fcplunge(i) = H./3./sqrt(3).*integral(int_part, 0, d(i));
    else
        Fcplunge(i) = Fcplungec;
    end
end

% 塑性阶段
Fcs = H./sqrt(3)./3.*As.*cot(pi/4-beta/2+Rake_angle_rad/2); % 公式第一项塑性阶段的切削力

% 公式第二项（不额外创建变量存储，直接最后加）
% 公式第三项（不额外创建变量存储，直接最后加）
% 公式最终结果
Fcp = Fcplunge+Fcs+H.*Ac./3+miuf.*Ff.*cos(alpha)-Ff.*sin(alpha);

% 切削能量
Edp = Fcp.*vc; % 塑性变形阶段所消耗的切削能量。直接用单位时间就好了，不用像周行那样还弄一个振动周期（因为我们不是用振动切削）(3-25)

% 塑性模式比切削能计算
Espd = Edp./Ac./vc./10^9; % 比切削能量(3-30)，单位用了GPa的感觉，总能量除以体积，因为我们没有使用振动周期，所以不用考虑振动频率f

%% 脆性阶段
k3 = 1; % 这个系数论文里也没说
v = Parameters.v; % 泊松比
cl = C2.*(1./tan(b)).^(5/12).*(E.^(3/4)./H./Kc./sqrt(1-v.^2)).^(1/2).*Ft.^(5/8); % 横向裂纹长度(3-36)
ch = C1.*(1./tan(b)).^(1/3).*E.^(1/2)./H.*Ft.^(1/2); % 横向裂纹深度(3-35)
cm = 0.12*(E/H)^(1/3).*cot(b).^(4/9)*(Ft./Kc).^(2/3); % 中间裂纹的长度(3-37)

Ab = 1/2*pi.*(ch).^2; % 脆性去除时垂直切削方向上材料去除的横截面积
Vb = Ab.*vc; % 脆性模式切削下材料去除率，周行论文为vb(3-31)
Ef = (2*pi.*cl+2.*cm).*vc.*ys; % 单位时间内横向裂纹和中间裂纹产生新表面所消耗的断裂能(3-38)

%% 相变所需能量
Aphase = Ab - Ac; % 脆性去除时硅的相变区域面积(整体面积Ah-塑性去除的面积Ac)，这个公式可能不太对，应该用仿真来获取刀以下多深的位置来计算相变面积
Mall = Aphase.*rou.*vc; % 相变区域的硅总质量，2.1e6为硅的密度
Ephase = Mall./Mphase.*Qc; % 相变能量.7e3应该是非晶硅晶化能，周行论文为Qc。28为硅的摩尔质量。(3-39)(3-40)

%% plastic zone的能量（切深以下塑性区域的塑性屈服能量）
sigmay = Parameters.sigmay; % 硅的屈服强度(Pa)
Vp = Aphase.*vc; % 脆性切削时切深以下的塑性区域体积
k4 = 1;
Ep = k4*Vp*sigmay; % 脆性切削时切深以下的塑性区域消耗的能量

%% 脆性阶段总能量
Eb = Ef+Edp+Ephase+Ep; % 脆性模式切削过程中单位时间消耗的总能量，Ef产生新表面所消耗的断裂能(3-38),Edp为塑性去除的能量(3-25),Ephase为相变能量(3-39),Ep为切深以下塑性区域的屈服能量,(3-41)

%% 脆性阶段比切削能量
Espb = Eb./Vb./10^9; % 在脆性模式中，相应比切削能量(3-42)

%% 找出Espd=Espb的点，即脆塑转变深度
critical_doc = 0; % 由于Ft的取值问题，可能导致Espb和Espd没有交点，所以人为设置0
for i = 1:numel(d)
    if abs(Espd(i)-Espb(i))<1
        critical_doc = d(i)*1E9;
        break;
    end
end

end