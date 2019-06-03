% 不要直接运行这个m文件，会跑很久很久很久很久……
% 需要哪部分复制到matlab运行就好

% 以下两行是Si和SiO2的折射率随波长的函数，可以先复制到matlab里运行
n_Si = @(lam) sqrt(11.6858+0.939816./(lam.*1e6).^2+0.000993358./((lam.*1e6).^2-1.22567));
n_SiO2 = @(lam) sqrt(1 + 0.6961663.*(lam.*1e6).^2./((lam.*1e6).^2-0.0684043.^2) + 0.4079426.*(lam.*1e6).^2./((lam.*1e6).^2-0.1162414.^2) + 0.8974794.*(lam.*1e6).^2./((lam.*1e6).^2-9.896161.^2));

% 将ModeSolver.m添加到matlab的工作目录

% ModeSolver是计算的核心，调用方式为
% ModeSolver(lambda, radius, n_core, n_cladding)
% 四个参数分别是波长(m)，半径(m)，芯层折射率，包层折射率
% 比如要计算lambda=1550nm, r=300nm, n1=1.44, n2=1, 则调用
% a = ModeSolver(1550e-9, 300e-9, 1.44, 1)
% 默认只计算HE11模
% 如果想计算其他模式，见下面的 可选参数 部分
% 结果在a.res里
% 比如想知道HE11模式的neff，调用a.res.HE11.neff()
% 想计算TE01模式的群速度，调用a.res.TE01.vg()

% 运行时会在工作目录生成一个名为tempFunctions的文件夹用于存放临时文件

% 可计算的物理量包括：
% 电场(r, theta)：E.r, E.x, E.y, E.z, E.phi, E.norm
% 磁场(r, theta)：H.r, H.x, H.y, H.z, H.phi, H.norm
% Sz(r, theta)
% 某一半径r以内的能量占总能量的比例Szint: Szint(r)
% 芯区的能量比例： eta (即Szint(radius))
% 群速度： vg
% 有效半径： Reff
% 色散： Dw

% 可选参数：
% HE11: 是否只计算HE11模，默认为true
% vvmax: 计算到HEvm或EHvm的阶数，比如想计算到HE31，则'vvmax', 3，想要计算到EH42，'vvmax', 4；默认为1；若该值大于1的话HE11自动设为false
% Dw: 是否计算色散(ture or false)，因为计算色散会比较慢，所以默认不计算，即'Dw', false

% 比如要计算lambda=1550nm, r=1000nm, n1=1.44, n2=1的情况下HE21模的色散，则调用
% a = ModeSolver(1550e-9, 1000e-9, 1.44, 1, 'Dw', true, 'vvmax', 2)
% a.res.HE21.Dw()

% 如果想要计算lam=1000~1500nm这种有自变量的情况，要用indeps参数

% 以下程序均对应2004OE-Single-mode guiding properties of subwavelength-diameter silica and silicon wire waveguides中的各图

% 例如fig 2，计算diameter=0~1000nm的情况
figure
lam = 633e-9;
Ds = linspace(0, 1000, 1001); % 单位是nm
n1 = n_SiO2(lam); % n_SiO2, n_Si输入的波长单位都应是m
a = ModeSolver(lam, @(D) D.*1e-9/2, n1, 1, 'indeps', {'D', Ds}, 'vvmax', 3);
% 最后的'indeps'参数指定了一个自变量'D', 值是Ds，即0~1000nm
% 前四个参数中radius受自变量D的影响，影响关系为 @(D) D.*1e-9/2 （单位从nm变成m再除以2），所以以匿名函数的形式写在这里
% 之后想要计算结果时要指定自变量的值，比如想计算'HE11'模式在D=500nm时的值，要写成a.res.HE11.neff(500)
% 也支持同时计算多个值，比如计算D＝[500, 501, 502]nm时的neff的值，可以写成a.res.HE11.neff([500, 501, 502])

modes = {'HE11', 'TE01', 'HE21', 'TM01', 'EH11', 'HE31', 'HE12', 'EH21'};
hold on
for ii=1:length(modes)
    plot(Ds, a.res.(modes{ii}).neff(Ds)*2*pi/(lam*1e6), 'DisplayName', modes{ii})
end
legend()

% fig 3
figure
lam = 1500e-9;
Ds = linspace(0, 800, 801);
n1 = n_Si(lam);
a = ModeSolver(lam, @(D) D.*1e-9/2, n1, 1, 'indeps', {'D', Ds}, 'vvmax', 3);
modes = {'HE11', 'TE01', 'HE21', 'TM01', 'EH11', 'HE31', 'HE12', 'EH21'};
hold on
for ii=1:length(modes)
    plot(Ds, a.res.(modes{ii}).neff(Ds)*2*pi/(lam*1e6), 'DisplayName', modes{ii})
end
legend()

% fig 5
figure
lam = 633e-9;
rs = [100, 200, 400, 800, 1600]/2*1e-9;
n1 = n_SiO2(lam);
Rs = linspace(0, 1.5, 1001)*1e-6;
as = {};
for ii=1:length(rs)
    r = rs(ii);
    as{ii} = ModeSolver(lam, r, n1, 1);
end
% fig5a
hold on
for ii=1:length(rs)
    r = rs(ii)
    a = as{ii}.res.HE11.E.r();
    % 这里想计算Er，先这么写一下返回a这个函数
    res = (Rs<r).*n1.^2.*a(Rs, pi/2)+(Rs>=r).*1.^2.*a(Rs, pi/2);
    % 这里再用a(Rs, pi/2)计算不同半径Rs, 同时theta=pi/2的Er的值
    plot(Rs*1e6, res./res(2), 'DisplayName', num2str(r*2*1e9));
end
legend()
% fig5b
hold on
for ii=1:length(rs)
    r = rs(ii)
    a = as{ii}.res.HE11.E.phi();
    res = double(a(Rs, pi/2));
    plot(Rs*1e6, res./res(2), 'DisplayName', num2str(r*2*1e9));
end
legend()
% fig5c
hold on
for ii=1:length(rs)
    r = rs(ii);
    a = as{ii}.res.HE11.E.z();
    res =1j*a(Rs, pi/2);
    plot(Rs*1e6, res, 'DisplayName', num2str(r*2*1e9));
end
legend()

% fig 6a
figure
lam = 633e-9;
r = 400e-9/2;
n1 = n_SiO2(lam);
a = ModeSolver(lam, r, n1, 1);
[R, theta] = meshgrid(0:0.01*r:3*r, 0:0.01*pi:2*pi);
X = R.*cos(theta);
Y = R.*sin(theta);
b = a.res.HE11.Sz();
mesh(X*1e6, Y*1e6, b(R, theta))

% fig 6b
figure
lam = 633e-9;
r = 200e-9/2;
n1 = n_SiO2(lam);
a = ModeSolver(lam, r, n1, 1);
[R, theta] = meshgrid(0:0.01*r:3*r, 0:0.01*pi:2*pi);
X = R.*cos(theta);
Y = R.*sin(theta);
b = a.res.HE11.Sz();
mesh(X*1e6, Y*1e6, b(R, theta))

% fig 7a
figure
lam = 633e-9;
Ds = linspace(0, 1.6, 1001);
n1 = n_SiO2(lam);
a = ModeSolver(lam, @(D) D.*1e-6/2, n1, 1, 'indeps', {'D', Ds});
plot(Ds, a.res.HE11.eta(Ds))

% fig 7b
figure
lam = 1500e-9;
Ds = linspace(0, 3, 1001);
n1 = n_SiO2(lam);
a = ModeSolver(lam, @(D) D.*1e-6/2, n1, 1, 'indeps', {'D', Ds});
plot(Ds, a.res.HE11.eta(Ds))

% fig 7c
figure
lam = 1500e-9;
Ds = linspace(0, 1.2, 1001);
n1 = n_Si(lam);
a = ModeSolver(lam, @(D) D.*1e-6/2, n1, 1, 'indeps', {'D', Ds});
plot(Ds, a.res.HE11.eta(Ds))

% fig 8a
figure
lam = 633e-9;
Ds = linspace(0.11, 0.8, 501);
n1 = n_SiO2(lam);
a = ModeSolver(lam, @(D) D.*1e-6/2, n1, 1, 'indeps', {'D', Ds});
semilogy(Ds, a.res.HE11.Reff(Ds)*1e6*2, '-', Ds, Ds, '--')

% fig 8b
figure
lam = 1500e-9;
Ds = linspace(0.3, 1.9, 501);
n1 = n_SiO2(lam);
a = ModeSolver(lam, @(D) D.*1e-6/2, n1, 1, 'indeps', {'D', Ds});
semilogy(Ds, a.res.HE11.Reff(Ds)*1e6*2, '-', Ds, Ds, '--')

% fig 8c
figure
lam = 1500e-9;
Ds = linspace(0.18, 0.8, 501);
n1 = n_Si(lam);
a = ModeSolver(lam, @(D) D.*1e-6/2, n1, 1, 'indeps', {'D', Ds});
semilogy(Ds, a.res.HE11.Reff(Ds)*1e6*2, '-', Ds, Ds, '--')

% fig 9a
figure
lam = 633e-9;
Ds1 = logspace(log10(0.15), 1, 201);
n1 = n_SiO2(lam);
a = ModeSolver(lam, @(D) D.*1e-6/2, n1, 1, 'indeps', {'D', Ds1});

lam = 1500e-9;
Ds2 = logspace(log10(0.3), 1, 201);
n1 = n_SiO2(lam);
b = ModeSolver(lam, @(D) D.*1e-6/2, n1, 1, 'indeps', {'D', Ds2});

semilogx(Ds1, a.res.HE11.vg(Ds1)./3e8, Ds2, b.res.HE11.vg(Ds2)./3e8)

% fig 9b
figure
lam = 1500e-9;
Ds = logspace(log10(0.15), 1, 401);
n1 = n_Si(lam);
a = ModeSolver(lam, @(D) D.*1e-6/2, n1, 1, 'indeps', {'D', Ds});
semilogx(Ds, a.res.HE11.vg(Ds)./3e8)

% fig 10a
figure
lams = linspace(0.3, 2.5, 1001);
rs = [200, 400, 600, 800, 1000, 1200]*1e-9/2;
n2 = 1;
as = {};
for ii=1:length(rs)
    r = rs(ii)
    as{ii} = ModeSolver(@(lam) lam.*1e-6, r, 0, 1, 'indeps', {'lam', lams}, 'n1lam', n_SiO2);
    % 由于波长为自变量，n1与波长有关，可以在参数n1lam中指定波长与折射率的关系，前面的n1可以随便填
end

hold on
for ii=1:length(rs)
    r = rs(ii);
    res = as{ii}.res.HE11.vg(lams);
    plot(lams, res/3e8, 'DisplayName', num2str(r*1e9*2))
end
xlim([0.3, 2.5])
legend

% fig 10b
figure
lams = linspace(1.2, 2.5, 1001);
rs = [200, 250, 300, 350, 400, 450]*1e-9/2;
n2 = 1;
as = {};
for ii=1:length(rs)
    r = rs(ii)
    as{ii} = ModeSolver(@(lam) lam*1e-6, r, 0, 1, 'indeps', {'lam', lams}, 'n1lam', n_Si);
end

hold on
for ii=1:length(rs)
    r = rs(ii);
    res = as{ii}.res.HE11.vg(lams);
    plot(lams, res/3e8, 'DisplayName', num2str(r*1e9*2))
end
xlim([1.2, 2.5])
legend

% fig 11a
figure
lam = 633e-9;
Ds1 = linspace(0.15, 1, 401);
n1 = n_SiO2(lam);
a = ModeSolver(lam, @(D) D*1e-6/2, n1, 1, 'indeps', {'D', Ds1}, 'Dw', true, 'n1lam', n_SiO2);

lam = 1500e-9;
Ds2 = linspace(0.3, 1, 401);
n1 = n_SiO2(lam);
b = ModeSolver(lam, @(D) D*1e-6/2, n1, 1, 'indeps', {'D', Ds2}, 'Dw', true, 'n1lam', n_SiO2);

plot(Ds1, a.res.HE11.Dw(Ds1)*1e6, Ds2, b.res.HE11.Dw(Ds2)*1e6)

% fig 11b
figure
lam = 1500e-9;
Ds = linspace(0.15, 0.5, 401);
n1 = n_Si(lam);
a = ModeSolver(lam, @(D) D*1e-6/2, n1, 1, 'indeps', {'D', Ds}, 'Dw', true, 'n1lam', n_Si);
plot(Ds1, a.res.HE11.Dw(Ds)*1e3)

% fig 12a
figure
lams = linspace(0.3, 2.5, 501);
rs = [200, 400, 600, 800, 1000, 1200]*1e-9/2;
n2 = 1;
as = {};
for ii=1:length(rs)
    r = rs(ii)
    as{ii} = ModeSolver(@(lam) lam*1e-6, r, 0, 1, 'indeps', {'lam', lams}, 'Dw', true, 'n1lam', n_SiO2);
end

hold on
for ii=1:length(rs)
    r = rs(ii);
    res = as{ii}.res.HE11.Dw(lams);
    plot(lams, res*1e6, 'DisplayName', num2str(r*1e9*2))
end
xlim([0.3, 2.5])
legend

% fig 12b
figure
lams = linspace(1.2, 2.5, 501);
rs = [200, 250, 300, 350, 400, 450]*1e-9/2;
n2 = 1;
as = {};
for ii=1:length(rs)
    r = rs(ii)
    as{ii} = ModeSolver(@(lam) lam*1e-6, r, 0, 1, 'indeps', {'lam', lams}, 'Dw', true, 'n1lam', n_Si);
end

hold on
for ii=1:length(rs)
    r = rs(ii);
    res = as{ii}.res.HE11.Dw(lams);
    plot(lams, res*1e3, 'DisplayName', num2str(r*1e9*2))
end
xlim([1.2, 2.5])
legend

