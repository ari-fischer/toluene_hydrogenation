function out = model_1(B);

global exp_table
kB = 1.380649E-23;
h = 6.62607015E-34;
R = 8.3144E-3; % kJ/mol


%assign parameter values to model variables
H1 = B(1);
S1 = B(2);
H_tol = B(3);
S_tol = B(4);

gam_tol = 6;
gam_TH5 = B(5);

%number of H-atoms added to MCHE to form the TS
x=2;
%number of H-atoms removed from toluene to form MASI
y=1;

c=6; %coordination number of the Pt 111 latice

eps = gam_TH5/gam_tol; %site ratio

%read the experimental data
Ts = exp_table(:,1);
rates_exp = exp_table(:,5);
p_H2 = exp_table(:,2);
p_tol = exp_table(:,3);
psi_exp = exp_table(:,6);

%edetermine coefficient values
A = kB.*Ts./h;
C1 = 1./(gam_tol.^eps).*A.*exp(-(H1-Ts.*S1)./R./Ts);
lamb = exp(-(H_tol-Ts.*S_tol)./R./Ts);

%solve polynomial to get H* site coverage (v_outs)
v_outs = [];
for i =1:length(p_tol)
F = @(v) v+gam_tol.*lamb(i).*p_tol(i)./(p_H2(i).^(y/2+gam_tol/2)).*v.^gam_tol-1;
v_out =  fsolve(F,.4);
v_outs = [v_outs,v_out];
end
v_outs=v_outs';

%psi equation 
r = c.*gam_TH5.*C1.*p_H2.^(x/2)./(p_H2).^(gam_TH5./2).*v_outs.^gam_TH5;

% calcuate normalized residuals as output
out = (r-psi_exp)./psi_exp;
