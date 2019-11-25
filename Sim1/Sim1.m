%Author: Kyle Felix
%Date: 10/25/19
%Class: ECE 416 Materials and Devices

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Problem 1a:  Plot energy vs. the occupation probability of electrons 
%in a material at T = 0, 150, 250, 350, and 450 K. Does the conductivity 
%increase/decrease with temperature? Why?
%This assumes the Fermi energy level is at the midpoint of Ec and Ev
clear;
k = 8.617e-5;           %Boltzman Constant
for ii = 1:5           %0-4 for loop 
    T = 100*(ii - 1);   %Sets the temp to 0,100, 200, 300, 400
    if ii >= 2
        T = T +50;      %Adds 50 to the temp -> 150,250,350,450 except for 0
    end    
    kT = k*T + 0.00000000000000000000000001;%Makes syntax easier for later on - Values checked
    dE(ii,1) = -5*kT;                       %Difference between E and EF
    for jj = 1:101                          %Each line on the graph is 100 points
        f(ii,jj) = 1/(1+exp(dE(ii,jj)/kT)); %Generic Fermi function
        dE(ii, jj+1) = dE(ii,jj) + 0.1*kT;  %Increasing E-Ef energy gap
    end 
end
dE = dE(:,1:jj);

close 
figure(1);
plot(dE',f'); grid;
xlabel('E-EF(eV)'); ylabel('f(E)');
text(.05, .2, 'T = 450K'); text(-.02, .05, 'T = 0K');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Problem 1b: Using the empirical fit equations* given for Si at 300 K, 
%plot electron and hole drift mobility (?n and ?p) as a function of doping 
%concentration (Nx) for Nx = 1014 to 1020 cm-3. Present on one graph.

%The following constants are used in the given equation on problem 1b
NDref = 6.95e18; NAref = 3.75e18;
unmin = 88; upmin = 54.3;
un0 = 1252; up0 = 407;
an = 1; ap = 1;

N = logspace(14,20);    %view the mobility ove the region 1e14 to 1e20
un = unmin + un0./(1+(N/NDref).^an);%Equation for donor doping
up = upmin + up0./(1+(N/NAref).^ap);%Equation for acceptor doping

%close 
figure(2);      
loglog(N, un, N, up); grid; %Disp;ay figure 3 with a log scale on both axis
axis([1.0e14 1.0e20 1.0e1 1.0e4]);  %Range of both axis
xlabel('NA or ND(cm-3)');
ylabel('Mobility(cm2/V-sec)');
text(1.0e15, 1500, 'Electrons');    %labels
text(1.0e15, 300, 'Holes');
text(1.0e18, 3000, 'Si at 300K');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Problem 1c: 

q = 1.6e-19;        %Constant magnitude of a charge of an electron  or hole
%Defining the mobility for each of the doping cases
un_0 = unmin + un0./(1+((1e10)/NDref).^an);
up_0 = upmin + up0./(1+((1e10)/NAref).^ap);
un_14 = unmin + un0./(1+((1e14+1e10)/NDref).^an);
up_14 = upmin + up0./(1+((1e14+1e10)/NAref).^ap);
un_16 = unmin + un0./(1+((1e16+1e10)/NDref).^an);
up_16 = upmin + up0./(1+((1e16+1e10)/NAref).^ap);
un_19 = unmin + un0./(1+((1e19+1e10)/NDref).^an);
up_19 = upmin + up0./(1+((1e19+1e10)/NAref).^ap);

%Solving for the resistivities for each of the doping cases
pa_0= 1/(q*1e10*up_0);
pa_14= 1/(q*(1e14+1e10)*up_14);
pa_16= 1/(q*(1e16+1e10)*up_16);
pa_19= 1/(q*(1e19+1e10)*up_19);
pd_0= 1/(q*1e10*un_0);
pd_14= 1/(q*(1e14+1e10)*un_14);
pd_16= 1/(q*(1e16+1e10)*un_16);
pd_19= 1/(q*(1e19+1e10)*un_19);

%Solving for the resistance given the resistivities calculated above
ra0= pa_0*75;
ra14= pa_14*75;
ra16= pa_16*75;
ra19= pa_19*75;
rd0= pd_0*75;
rd14=pd_14*75;
rd16=pd_16*75;
rd19=pd_19*75;

%This is the range and step size of voltages for the plot
V = 0:0.001:1;

%Calculating the current for each of the doping cases
ia0= V/ra0;
ia14= V/ra14;
ia16=V/ra16;
ia19=V/ra19;
id0= V/rd0;
id14=V/rd14;
id16=V/rd16;
id19=V/rd19;

%close
figure(3);
%Display a plot with log scaled x and y axis with all of the cases
loglog(V, ia0, V, ia14, V, ia16, V, ia19, V, id0, V, id14, V, id16, V, id19); 
grid on;
xlabel('Voltage (V)');
ylabel('Current (I)');

%tabulate the resistivities in the Command Window
table(pa_0,pa_14,pa_16,pa_19,pd_0,pd_14,pd_16,pd_19)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




