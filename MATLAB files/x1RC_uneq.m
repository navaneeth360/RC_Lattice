%Finding the objective function
%When we have a 2 sector lattice would unequal capacitances help

n = 1;      %No. of stages
R = 1;    %Total resistance
C = 1;  %Each capacitor
T_d = 2*R/2*C;


f_min = 0.0001; 
f_max = 1.1/8;

num_points_per_decade = 250;
num_decades = log10(f_max) - log10(f_min);
num_points = num_points_per_decade * num_decades;
frequencies = linspace(f_min, f_max, 1000);

w = 2*pi*frequencies;
s = 1i*w;

r_mid = R/4;
func2 = @(x) sinc(x).*exp(-1j*2*pi.*x);
func1 = @(r1, r2, r3, x, C) 2./( (r2 + 2*r3) + (2*r1 + r2)*(4*r3 + r2 + (1./(1j*2*pi*C.*x)))./((1./(1j*2*pi*C.*x) - r2)));
func3 = @(x,C) (1 - 1j*2*pi*0.5*C.*x)/(1 + 1j*2*pi*0.5*C.*x);

%objective = @(u) integral(@(x) 1000*(abs(func1(u(1),u(2),x,u(3)) - func2(x)).^2),0,1.1/8);
%objective2 = @(u) integral(@(x) (abs(func1(0.25,0.5,x,u) - func2(x)).^2),0,1.1/8);
%objective2 = @(u) integral(@(x) (abs(func1(0.25,0.5,x,u) - func2(x)).^2),0,1.1/8);

objective3 = @(u) integral(@(x) (abs(func1(u(1), u(2), u(3), x, u(4)) - func2(x)).^2),0,1.1/8);
objective4 = @(u) max(abs(func1(u(1), u(2), u(3), frequencies, u(4)) - func2(frequencies)));

options = optimoptions('fminunc', 'Display', 'iter', 'Algorithm', 'quasi-newton');
r1_0 = 0.25;
r2_0 = 0.5;
r3_0 = 0.25;
c_0 = 1;
%[k_opt, fval] = fminunc(objective, k0, options);

Aeq = [1, 1, 1, 0];
beq = 1;
lb = [0.00001 0.00001 0.00001 0];
ub = [];
x0 = [r1_0 r2_0 r3_0 c_0];
x00 = [0 0 0 0]
x1 = [0.2731    0.4539    0.2731    1.0492];
x2 = [0.2719    0.4561    0.2719    1.0382];

opti3 = fmincon(objective3, x1, [],[],Aeq,beq,lb,ub,[],options);
opti5 = fmincon(objective4, x1, [],[],Aeq,beq,lb,ub,[],options);
opti4 = fmincon(objective4, x1, [],[],Aeq,beq,lb,ub);


disp(opti3)
disp(opti5)

disp(objective3(opti3))
disp(objective3(opti4))
%disp(objective4([0.25 0.5 0.25 1.0042]))

for i = 1:1000
    arr41(i) = (angle(func1(opti3(1),opti3(2),opti3(3),frequencies(i),opti3(4)))*180/pi);
    arr421(i) = abs(func1(opti3(1),opti3(2),opti3(3),frequencies(i),opti3(4)) - func2(frequencies(i)));
    arr431(i) = abs(func1(opti5(1),opti5(2),opti5(3),frequencies(i),opti5(4)) - func2(frequencies(i)));
    arrx(i) = angle(func1(x2(1),x2(2),x2(3),frequencies(i),x2(4)))*180/pi;
end

hold on
%plot(frequencies, arr421, LineWidth=3)
plot(frequencies, arr431, LineWidth=3)
%plot(frequencies, arr4200, LineWidth=3)
plot(frequencies, arr4201, LineWidth=3)
%plot(frequencies, arr430, LineWidth=3)
plot(frequencies, arr43, LineWidth=3)
hold off
%legend("1x Optimized 2 norm","1x Optimized infinity norm","2x Optimized 2 norm","2x Optimized infinity norm","3x Optimized 2 norm","3x Optimized infinity norm",fontsize=20)
legend("1x Optimized infinity norm","2x Optimized infinity norm","3x Optimized infinity norm",fontsize=20)
xlabel('Frequency (in Hz)',FontSize=20);
ylabel('Magnitude of Difference',FontSize=20);
%ylabel('Phase (in degrees)',FontSize=20);
title('Magnitude of Difference for optimum',FontSize=20)
%title('Phase Response',FontSize=20)
%title('Difference in magnitude for optimum',FontSize=20)
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
grid on
writematrix(arr41, 'uneq12.dat')
