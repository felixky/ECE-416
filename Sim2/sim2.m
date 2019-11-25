%Author: Kyle Felix
%Date: 11/20/19
%Description: Sim2 of GaAs

clc
clear all
close all
%Constants
ni = 2.1e6;
Na = [1e19, 5e17, 1e16, 5e15];
Nd = [1e15, 1e16, 5e17, 1e18];
T = 300;
tau = 15e-6;
taud = 100e-6;
Kgaas = 12.9;
E0 = 8.85e-14;
H = 0.00003;   %3um
L = 0.0002;   %2um
WL = 0.00015;  %1.5um
A = L*WL;
Va = -5;
m = 0.5;
k = 8.617e-5;
q = 1.6e-19;
%Mobility, diffusivity, and lifetime
un = 500 + ((8900)./(1 + (Nd./(6e16))).^0.394);
up = 20+471.5./(1+Na./(1.48E17)).^(0.38);
Dp = (0.0259).*up;
Dn = (0.0259).*un;
Ln = (tau.*Dn).^(1/2);
Lp = (tau.*Dp).^(1/2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Problem 1a:Assuming ideal diode operation, calculate and tabulate
%(make a table of) the following values for each doping profile listed above:
% I0: Reverse saturation current [in A]
% Cj0: zero-bias junction capacitance [in F]
% Cj, -5V: junction capacitance at -5V bias (5V reverse bias) [in F]
% Vbi: built-in voltage [in V]
Vbi = (0.0259).*log((Na.*Nd)./(ni^2));
p = table(Vbi(1),Vbi(2),Vbi(3),Vbi(4));
p.Properties.VariableNames = {'Vbi_i','Vbi_ii','Vbi_iii','Vbi_iv'};
p

W = (( ((2*Kgaas*E0)/q) .* ((Na+Nd)./(Na.*Nd)) .*(Vbi-Va)) .^(1/2));

Cj0 = (Kgaas*E0*A)./(( ((2*Kgaas*E0)/q) .* ((Na+Nd)./(Na.*Nd)) .*(Vbi)) .^(1/2));
p = table(Cj0(1),Cj0(2),Cj0(3),Cj0(4));
p.Properties.VariableNames = {'Cj0_i','Cj0_ii','Cj0_iii','Cj0_iv'};
p

Cj = (Kgaas*E0*A)./W;
p = table(Cj(1),Cj(2),Cj(3),Cj(4));
p.Properties.VariableNames = {'Cj_i','Cj_ii','Cj_iii','Cj_iv'};
p

I0 = (q*A).*( ( (Dn.*(ni^2)) ./ (Ln.*Na) ) + ( (Dp.*(ni^2)) ./ (Lp.*Nd) ) );
p = table(I0(1),I0(2),I0(3),I0(4));
p.Properties.VariableNames = {'I0_i','I0_ii','I0_iii','I0_iv'};
p
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Problem 1b:Using these values, create log(|I|)-V curves from -1.5 to 1.5 V
%for each doping profile. Present all four plots on one graph, labeling 
%everything clearly and including positive and negative voltages. Again, 
%assume ideal diode operation and note that the absolute value is needed 
%for a log plot.

V = linspace(-1.5,1.5);
Irev = zeros(length(I0), length(V));
 for x = 1:length(I0)         %Calculating the (W)idth of the depletion region
    Irev(x,:) = I0(x).*(exp((V)./(0.0259))-1);
 end

figure(1);
semilogy(V,abs(Irev)); grid;
xlabel('Voltage (V)'); ylabel('Current (A)');
title('I-V for different Na and Nd');
legend('Na=1e19 Nd=1e15','Na=5e17 Nd=1e16','Na=1e16 Nd=5e17','Na=5e15 Nd=1e18');
%text(-1, 1e-25, 'Na=1e19 Nd=1e15');text(.05, .2, 'T = 450K');text(.05, .2, 'T = 450K');text(.05, .2, 'T = 450K');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Problem 1d:Compute the I-V curve for doping profile iii) with a minority 
%carrier lifetime of 100?s. Plot both I-V curves (15?s and 100?s) on a 
%graph. Explain what differences there are and why.
Lnd = (taud.*Dn).^(1/2);
Lpd = (taud.*Dp).^(1/2);
I0d = (q*A).*( ( (Dn.*(ni^2)) ./ (Lnd.*Na) ) + ( (Dp.*(ni^2)) ./ (Lpd.*Nd) ) );

Irevd = zeros(length(I0d), length(V));
 for x = 1:length(I0d)         %Calculating the (W)idth of the depletion region
    Irevd(x,:) = I0d(x).*(exp((V)./(0.0259))-1);
 end

figure(2);
semilogy(V,abs(Irev(3,:)), V,abs(Irevd(3,:))); grid;
xlabel('Voltage (V)'); ylabel('Current (A)');
title('Differing Minority Lifetimes');
legend('Tau = 15uS Na=1e16 Nd=5e17', 'Tau = 100uS Na=1e16 Nd=5e17');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
