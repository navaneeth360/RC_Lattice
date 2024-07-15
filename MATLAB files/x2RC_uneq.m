%Finding the objective function
%When we have a 2 sector lattice would unequal capacitances help

n = 1;      %No. of stages
R = 1;    %Total resistance
C = 1;  %Each capacitor
T_d = 2*R/2*C;

R6 = @(r1,r2,r3,r4,r5) 1 - r1 - r2 - r3 - r4 - r5;

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

%V2/I6, I6 is I_out
tf1 = @(r4, r5, r6, c2, x) (r5 + 2*r6) + (2*r4 + r5)*(4*r6 + r5 + (1./(1j*2*pi*c2.*x)))./(1./(1j*2*pi*c2.*x) - r5);
% (2r3 + r2)i3/i6
tf2 = @(r2, r3, r5, r6, c2, x) (2*r3 + r2)*(4*r6+ r5 + (1./(1j*2*pi*c2.*x)))./((1./(1j*2*pi*c2.*x)) - r5);
% (2r1 + r2)i1/i6
tf3 = @(r1, r2, r3, c, r5, r6, c2, x) (2*r1 + r2)*(4*r3 + r2 + (1./(1j*2*pi*c.*x)))./((1./(1j*2*pi*c.*x)) - r2).*(4*r6+ r5 + (1./(1j*2*pi*c2.*x)))./((1./(1j*2*pi*c2.*x)) - r5);
%Extra term from i1 and i3 relationship
tf4 = @(r1, r2, r4, r5, r6, c, c2, x) 2*(2*r1 + r2)./((1./(1j*2*pi*c.*x)) - r2).*tf1(r4, r5, r6, c2, x);


func1 = @(r1, r2, r3, r4, r5, r6, c, c2, x) 2./(tf1(r4, r5, r6, c2, x) + tf2(r2, r3, r5, r6, c2, x) + tf3(r1, r2, r3, c, r5, r6, c2, x) + tf4(r1, r2, r4, r5, r6, c, c2, x));

objective = @(u) integral(@(x) (abs(func1(u(1),u(2),u(3),u(4),u(5),u(6),u(7),u(8),x) - func2(x)).^2),0,1.1/8);

objective2 = @(u) max(abs(func1(u(1),u(2),u(3),u(4),u(5),u(6),u(7),u(8),frequencies) - func2(frequencies)));

options = optimoptions('fminunc', 'Display', 'iter', 'Algorithm', 'quasi-newton');
r1_0 = 1;
r2_0 = 1;
r3_0 = 1;
r4_0 = 1;
r5_0 = 1;
r6_0 = 1;
c1_0 = 1;
c2_0 = 1;
x1 = [ 0.0429    0.2279    0.2292    0.2292    0.2279    0.0429    0.9595    0.9595];
% 8.9243e-07 - final
x2 = [0.5 0.5 0.1 0.1 0.1 0.5 1 1];
x3 = [0.0427    0.2277    0.2296    0.2296    0.2277    0.0427    0.9578    0.9578];
%8.9324e-07

x0 = [r1_0 r2_0 r3_0 r4_0 r5_0 r6_0 c1_0 c2_0];
A = [1, 1, 1, 1, 1, 1, 0, 0];
b = 1;
lb = [0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 ];
ub = [1 1 1 1 1 1 2 2];

opti1 = fmincon(objective, x3, [],[],A,b,lb,ub,[],options);
opti3 = fmincon(objective2, x3, [],[],A,b,lb,ub,[],options);
opti2 = fmincon(objective2, x3, [],[],A,b,lb,ub);
%opti3 = fmincon(objective, x0, [],[],A,b,lb,ub,options2);


for i = 1:1000
    arr41(i) = angle(func1(opti1(1),opti1(2),opti1(3),opti1(4),opti1(5),opti1(6),opti1(7),opti1(8),frequencies(i)))*180/pi;
    arr4200(i) = (abs( func1(opti1(1),opti1(2),opti1(3),opti1(4),opti1(5),opti1(6),opti1(7),opti1(8),frequencies(i)) - func2(frequencies(i)) ));
    arr4201(i) = (abs( func1(opti3(1),opti3(2),opti3(3),opti3(4),opti3(5),opti3(6),opti3(7),opti3(8),frequencies(i)) - func2(frequencies(i)) ));
    arrx(i) = angle(func1(x1(1),x1(2),x1(3),x1(4),x1(5),x1(6),x1(7),x1(8),frequencies(i)))*180/pi;
end

hold on
plot(frequencies, arr4200, LineWidth=3)
hold off
legend("Optimized",fontsize=20)
xlabel('Frequency (in Hz)',FontSize=20);
%ylabel('Magnitude of Difference',FontSize=20);
ylabel('Phase (in degrees)',FontSize=20);
%title('Difference in magnitude for optimum',FontSize=20)
title('Phase Response',FontSize=20)
%title('Difference in magnitude for optimum',FontSize=20)
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
grid on


%opti2 = fmincon(objective, x0, A,b,[],[],lb,ub);
disp(opti1)
disp(opti3)


disp(objective(opti1))
disp(objective(opti3))

disp(objective2([1/8 1/4 1/8 1/8 1/4 1/8 1.04 1.04]))
writematrix(arr41, 'uneq21.dat')
