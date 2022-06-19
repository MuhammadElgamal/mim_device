%% Full Integration of Code
clear all;
clc;
close all;

device_type = 'step';   % can be changed to resonant
voltage_limit = 0.5;
num_of_points = 25;
plotfigures=false;

if strcmp(device_type,'step')
    display('Step Tunneling Diode')
    aff=[4.29 3.88];
elseif strcmp(device_type,'resonant')
    display('Resonant Tunneling Diode')
    aff=[4.26 3.84];
end
% Source Metals Work Function
% for step=4.558
Left_metal_WF = 4.558;
% Drain Metal WF for NbO3
% for step=4.55
Right_metal_WF = 4.55;
% Ef left metal val
Ef_Left = 10;
% Channel Dielectrics
e1=[ 25 20 ];
% Oxide Permitivity
e2=[60 60 60];
% Channel Electrons Affinities for Step MIIM at 2nm


%% ------------<< Widths >>------------- %%
% ---- Nb2O5(4.23/25) - Ta2O5(3.83/20) ------%
%for Step (2nm)
w1=[ 1 1 ];
w2=[0.5 1 0.5]; % Neglect for diodes

%% ------------<< Heights >>------------ %%
h1=2;
h2=[1,3,1];

%-----------------------------------------------------------%
% Device Building with any initial potential distribution
device=mim_device('diode', {w1,w2}, {h1,h2}, {e1, e2},aff,300);
setMetalsSpecs(device,Ef_Left,Left_metal_WF,Right_metal_WF);

%% -------------------------------------------------------%
% Device Building with any initial potential distribution
device=mim_device('diode', {w1,w2}, {h1,h2}, {e1, e2},aff,300);
setMetalsSpecs(device,Ef_Left,Left_metal_WF,Right_metal_WF);
%% Calcualte IV Characteristics
device.VoltageLimit=voltage_limit;
device.VoltagePoints=num_of_points;
variable_to_fix='vgs';  fixed_value=0; plotfigures=false; slice_count=2;
eval("characterise(device, variable_to_fix, fixed_value,plotfigures, slice_count)");

%% Plotting Current and Other Charcteristics
% Run this part to show calibrartion
load('groover_data.mat')
% [20 100 100] are smoothing factors for current, first derivative of
% current, and second derivative of current
device.extract_features('TLM', [20 100 100]);
figure;
semilogy(device.volt, abs(device.curr));
title('IV Charctersistics')
hold on;
if strcmp(device_type,'step')
    plot(Groover.step.x, Groover.step.y, '--k')
    save('Step_current.fig')
elseif strcmp(device_type,'resonant')
    semilogy(Groover.resonant.x, abs(Groover.resonant.y), '--k')
    save('Resonant_current.fig')
end
legend('Ours', 'Grover');
figure;
semilogy(device.volt, device.resistivity*1.5e10);
hold on;
plot(groover.volt, log10(groover.resistivity*1.5e10), '--k');
title('Resistivity')

figure;
plot(device.volt, device.responsivity);
hold on;
plot(groover.volt, groover.responsivity, '--k');
title('Responsivity')
    
figure;
semilogy(device.volt, device.nonlinearity);
hold on;
plot(groover.volt, groover.nonlinearity, '--k');
title('Absolute Nonlinearity')

figure;
semilogy(device.assymetry.x, device.assymetry.y);
hold on;
plot(groover.assymetry.x, groover.assymetry.y);
title('Absolute Assymetry')
