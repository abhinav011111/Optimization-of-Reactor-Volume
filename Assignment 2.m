Ca0 = [2 5 6 6 11 14 16 24];
Caf = [0.5 3 1 2 6 10 8 4];
t = [30 1 50 8 4 20 20 4];
Ca_in = 10;
Ca_out = 1;
Xa_f = 0.9;

Xa = (Ca0 - Caf) ./ Ca0;
r_a = (Ca0 - Caf) ./ t;
rate_inv = 1 ./ r_a;

f = spline(Caf, rate_inv);
x = linspace(min(Caf), max(Caf), 1000);
y = ppval(f, x);

plot(x,y,"LineWidth",2);
hold on;
scatter(Caf,rate_inv);
xlabel('Caf');
ylabel('-1/ra');
title('Spline Interpolation of -1/ra vs Ca');

Volumetric_rate= 0.1;
% Volume of PFR
Volume_PFR = Volumetric_rate * integral(@(x) ppval(f, x), Ca_out, Ca_in);

%Volume of CSTR
Volume_CSTR = Volumetric_rate * (Ca_in - Ca_out) * ppval(f, Ca_out);

%Volume of PFR and MFR system
Ca_min = fminbnd(@(x) ppval(f, x), Ca_out, Ca_in);
rate_min = ppval(f, Ca_min);
Volume_PFR_MFR = Volumetric_rate * (integral(@(x) ppval(f, x), Ca_out, Ca_min) + (Ca_in - Ca_min) * rate_min);

%2 MFRs 
ymaxi=0;
xmaxi=0;
for x = 1 : 0.01 : 10
    area = (Ca_in - x)*(ppval(f,Ca_out) - ppval(f,x));
    if area > ymaxi
        ymaxi = area;
        xmaxi = x;
    end
end

Area1 = (xmaxi - Ca_out)*(ppval(f,Ca_out));
Area2 = (Ca_in - xmaxi)*(ppval(f,xmaxi));
Volume_2_MFRs = Volumetric_rate*(Area1 + Area2);

%PFR with Recycle
for x = 1.1 : 0.01 : 10
    area_pfr = integral(@(x) ppval(f, x), Ca_out,x);
    area_mfr = (x - Ca_out)*ppval(f,x);

    if abs(area_pfr - area_mfr) < 0.01
        Cai = x;
        break;
    end
end

R = (Ca_in - Cai)/(Cai-Ca_out);
Area_R = (Ca_in - Ca_out)*ppval(f,Cai);
Volume_PFR_Recycle = Area_R*Volumetric_rate;


fprintf('Single PFR                       : %f\n', Volume_PFR);
fprintf('Single CSTR                      : %f\n', Volume_CSTR);
fprintf('Two Stirred Tanks                : %f\n', Volume_2_MFRs);
fprintf('Combination of a PFR and a MFR   : %f\n', Volume_PFR_MFR);
fprintf('PFR with recycle                 : %f\n', Volume_PFR_Recycle);
