%% System identification script to fit a stat-space model to frequency response data of the circuit
% Extraction of state-space model of transformer
clear all
close all

%% Import frequency response from LT Spice export (created with AC sweep)

% Import waveforms form LTSpice
% Set up the Import Options
opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["frequency", "Vout", "Vw1", "I2", "I1", "Iw1", "Iw2"];
opts.VariableTypes = ["double", "string", "string", "string", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import sweep
dataTable = readtable("Transformer_FrequencyResponse.txt", opts);

f = dataTable.frequency;
df = f(2)-f(1);
tempVout = dataTable.Vout;
tempI1 = dataTable.I1;
tempI2 = dataTable.I2;
tempVw1 = dataTable.Vw1;
tempIw1 = dataTable.Iw1;
tempIw2 = dataTable.Iw2;
% convert strings into complex numbers 
for i=1:length(tempVout)
    Vout(i) = str2number(tempVout(i));
    Iw1(i) = str2number(tempIw1(i));
    Iw2(i) = str2number(tempIw2(i));
    Vw1(i) = str2number(tempVw1(i));
    I1(i) = str2number(tempI1(i));
    I2(i) = str2number(tempI2(i));
end

Vout = transpose(Vout);
Iw1 = transpose(Iw1);
Iw2 = transpose(Iw2);
Vw1 = transpose(Vw1);
I1 = transpose(I1);
I2 = transpose(I2);

Ts = 0;

% Create frequency response data objects
data_Vout = idfrd(Vout, f, Ts, 'FrequencyUnit','Hz');
data_Iw1 = idfrd(Iw1, f, Ts, 'FrequencyUnit','Hz');
data_Iw2 = idfrd(Iw2, f, Ts, 'FrequencyUnit','Hz');
data_Vw1 = idfrd(Vw1, f, Ts, 'FrequencyUnit','Hz');
data_I1 = idfrd(I1, f, Ts, 'FrequencyUnit','Hz');
data_I2 = idfrd(I2, f, Ts, 'FrequencyUnit','Hz');


% Combine into 2 data sets
% terminal signals
ResponseData_combined = [];
ResponseData_combined(1,1,:) = Vout;
ResponseData_combined(2,1,:) = I1;
ResponseData_combined(3,1,:) = I2;

data_combined = idfrd(ResponseData_combined, f, Ts, 'FrequencyUnit','Hz'); 

%winding signals
ResponseDataWinding_combined = [];
ResponseDataWinding_combined(1,1,:) = Vw1;
ResponseDataWinding_combined(2,1,:) = Iw1;
ResponseDataWinding_combined(3,1,:) = Iw2;

dataWinding_combined = idfrd(ResponseDataWinding_combined, f, Ts, 'FrequencyUnit','Hz'); 


%% Sys ID
Options = ssestOptions;                           
Options.Display = 'on'; 

nx = 1:10;

ssCombined = ssest(data_combined, nx);                                      % n = 4
figure()
compare(data_combined, ssCombined)

ssWindingCombined = ssest(dataWinding_combined, nx);                        % n = 4
figure()
compare(dataWinding_combined, ssWindingCombined)

%% Verification with ground truth (GT), exported from LT Spice (e.g. step response)
% Set up the Import Options
opts = delimitedTextImportOptions("NumVariables", 8);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["time", "Vout", "Vin", "Vw1", "Iw1", "Iw2", "I1", "I2"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% import data for verification:
dataTable = readtable("Transformer_VerificationData.txt", opts);
t_sim = dataTable.time;
Vout_GT = dataTable.Vout;
Vin_GT = dataTable.Vin;
Vw1_GT = dataTable.Vw1;
Iw1_GT = dataTable.Iw1;
Iw2_GT = dataTable.Iw2;
I1_GT = dataTable.I1 * (-1);
I2_GT = dataTable.I2 * (-1);

%% Sim and Plot
V_in_test = Vin_GT;

combined_sim = lsim(ssCombined,V_in_test,t_sim);
Windingcombined_sim = lsim(ssWindingCombined,V_in_test,t_sim);

figure()

subplot(3,2,1)
hold on
grid on
plot(t_sim, combined_sim(:,1))
plot(t_sim, Vout_GT)
legend('SysID', 'GT')
title('v_{2}')
xlabel('Time in (s)')
ylabel('Voltage in (V)')

subplot(3,2,3)
hold on
grid on
plot(t_sim, combined_sim(:,2))
plot(t_sim, I2_GT)
legend('SysID', 'GT')
title('i_1')
xlabel('Time in (s)')
ylabel('Current in (A)')

subplot(3,2,5)
hold on
grid on
plot(t_sim, combined_sim(:,3))
plot(t_sim, I1_GT)
legend('SysID', 'GT')
title('i_2')
xlabel('Time in (s)')
ylabel('Current in (A)')

subplot(3,2,2)
hold on
grid on
plot(t_sim, Windingcombined_sim(:,2))
plot(t_sim, Iw1_GT)
legend('SysID', 'GT')
title('i_{w1}')
xlabel('Time in (s)')
ylabel('Current in (A)')

subplot(3,2,4)
hold on
grid on
plot(t_sim, Windingcombined_sim(:,3))
plot(t_sim, Iw2_GT)
legend('SysID', 'GT')
title('i_{w2}')
xlabel('Time in (s)')
ylabel('Current in (A)')

subplot(3,2,6)
hold on
grid on
plot(t_sim, Windingcombined_sim(:,1))
plot(t_sim, Vw1_GT)
legend('SysID', 'GT')
title('v_{w1}')
xlabel('Time in (s)')
ylabel('Voltage in (V)')