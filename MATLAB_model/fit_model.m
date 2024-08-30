clc
clear
close all

global exp_table
%import kinetic data
file = 'data_psi_HT_full.csv'
exp_table = table2array(readtable(file))

%initial parameter values
B = [-49.9375167123973	-0.269657887644013	-38.5104700556806	-0.0392038886705180, 4];

%parameter ranges
lb = [-300, -0.5...
    -200, -.3,1];
ub = [100, 0.1,...
    100,.1,10];

%run regression analysis
options = optimoptions('lsqnonlin','Display','iter')
[beta_fit,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(@model_1,B,lb,ub, options)%,xdata,ydata, [0], [10], options)

%calculate confidence intervals
ci = nlparci(beta_fit,residual,'jacobian',jacobian)

%get the final normalized residual values
outs = model_1(beta_fit)

%parity plot
figure(1)
hold on
plot(exp_table(:,6),outs.*exp_table(:,6)+exp_table(:,6),'x')
plot ((1:80000),(1:80000))
hold off

%Mean absolute percentage errors
mean(abs(outs.*exp_table(:,6)))./mean(exp_table(:,6))
%export results
save('output_variables')



