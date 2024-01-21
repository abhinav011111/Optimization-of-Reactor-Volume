k1 =1;
k2=2;
k3=3;
syms Ca 
f = k1 * Ca ./ (1 + k2 * Ca + k3*Ca^2);

% calculate derivative
df = diff(f,Ca);
Ca1=2;
df_2 = subs(df,Ca,Ca1);
df_2 = double(df_2);

% calculate area
area = int(f,0,4);
area = double(area);

%max of function
df0 = df==0;
max_Ca = solve(df0,Ca);
max_Ca = max(max_Ca);
max_Ca = double(max_Ca);
max_f = subs(f,Ca,max_Ca);
max_f = double(max_f);

%line between two points
Ca2 = 0.4;
Ca3 = 2.5;
points_range = [Ca2 Ca3];
f2 = subs(f,Ca,Ca2);
f3 = subs(f,Ca,Ca3);
slope = (f3-f2)/(Ca3-Ca2);
slope = double(slope);

%point where the tangent is the same as the given slope
df1 = df==slope;
assume(Ca,'real');
point = vpasolve(df1,Ca,[0,4]);
point = double(point);
f_point = subs(f,Ca,point);
f_point = double(f_point);
y = slope*(Ca-point) + f_point;
%-ra VS Ca
fplot(f,[0 4],"LineWidth",2,"Color","#A2142F");
hold on
%line between two points
plot(points_range, [f2 f3],"LineStyle","-","LineWidth",2,"Color",'m');
hold on
plot(max_Ca,max_f,'ob:', 'MarkerFaceColor','k');
hold on
plot(point,f_point,'ob:', 'MarkerFaceColor','r');
hold on
%tangent
fplot(y,[0 10],"LineWidth",2,"Color","#80B3FF");
hold off
legend('Rate expression','Line between Ca = 1 and Ca = 4','max_rate', 'Tangent at point','Derivative at Ca = 2')
xlabel('Ca')
ylabel('-ra')
title('Plot of -ra VS Ca')