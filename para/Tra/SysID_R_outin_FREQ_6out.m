%% Linear System identification
% Extraction of state-space model of transformer
% Comparison of Rc inside and outside in equivalent circuit for testing
clear all
close all
%% Frequency response via AC sweep in LT Spice

%% Import waveforms form LTSpice
%% Set up the Import Options
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


%% Import sweep (dec)
%dataTable = readtable("C:\Users\q634145\OneDrive - BMW Group\04_Praxiseinsätze\ES-663\04_LTSpice\TransformerR_Out_FREQ_dec", opts);
%dataTable = readtable("C:\Users\q634145\OneDrive - BMW Group\04_Praxiseinsätze\ES-663\04_LTSpice\TransformerR_In_R1e3_dec", opts);
% Capacitor Test
%dataTable = readtable("C:\Users\q634145\OneDrive - BMW Group\04_Praxiseinsätze\ES-663\04_LTSpice\System Identification\FREQ_Cap\Transformer_3Cap_Sweep_6outputs_150pF", opts);
% Measured Transformer
dataTable = readtable("C:\Users\q634145\OneDrive - BMW Group\04_Praxiseinsätze\ES-663\04_LTSpice\System Identification\MeasuredTransformer\PlanarTransformer_Sweep_Final.txt", opts);

f = dataTable.frequency;
df = f(2)-f(1);
tempVout = dataTable.Vout;
tempI1 = dataTable.I1;
tempI2 = dataTable.I2;
tempVw1 = dataTable.Vw1;
tempIw1 = dataTable.Iw1;
tempIw2 = dataTable.Iw2;
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

data_Vout = idfrd(Vout, f, Ts, 'FrequencyUnit','Hz');
data_Iw1 = idfrd(Iw1, f, Ts, 'FrequencyUnit','Hz');
data_Iw2 = idfrd(Iw2, f, Ts, 'FrequencyUnit','Hz');
data_Vw1 = idfrd(Vw1, f, Ts, 'FrequencyUnit','Hz');
data_I1 = idfrd(I1, f, Ts, 'FrequencyUnit','Hz');
data_I2 = idfrd(I2, f, Ts, 'FrequencyUnit','Hz');


ResponseData_combined = [];
ResponseData_combined(1,1,:) = Vout;
ResponseData_combined(2,1,:) = I1;
ResponseData_combined(3,1,:) = I2;
%ResponseData_combined(3,1,:) = Vw1;
%ResponseData_combined(4,1,:) = Iw1;
%ResponseData_combined(5,1,:) = Iw2;

data_combined = idfrd(ResponseData_combined, f, Ts, 'FrequencyUnit','Hz'); 

ResponseDataWinding_combined = [];
ResponseDataWinding_combined(1,1,:) = Vw1;
ResponseDataWinding_combined(2,1,:) = Iw1;
ResponseDataWinding_combined(3,1,:) = Iw2;

dataWinding_combined = idfrd(ResponseDataWinding_combined, f, Ts, 'FrequencyUnit','Hz'); 
% figure()
% bode(data_Vout)
% figure()
% bode(data_Iw1)
% figure()
% bode(data_Iw2)


%% Sys ID
Options = ssestOptions;                           
Options.Display = 'on'; 

nx = 1:10;
% ssVo = ssest(data_Vout, nx);    %4
% figure()
% compare(data_Vout, ssVo)
% ssIw1 = ssest(data_Iw1, nx);    %5
% figure()
% compare(data_Iw1, ssIw1)
% ssIw2 = ssest(data_Iw2, nx);    %4
% figure()
% compare(data_Iw2, ssIw2)

ssCombined = ssest(data_combined, nx);                                      % 200pF: nx = 6, 150pF: nx = 6, measuredPara: 6, Final: 4
figure()
compare(data_combined, ssCombined)
%%
ssWindingCombined = ssest(dataWinding_combined, nx);                        % 200pF: nx = 6, 150pF: nx = 4, measuredPara: 4, Final: 4
figure()
compare(dataWinding_combined, ssWindingCombined)

%% Verification
%% Set up the Import Options
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

% Data for verification:

% Capacitor Test
%dataTable = readtable("C:\Users\q634145\OneDrive - BMW Group\04_Praxiseinsätze\ES-663\04_LTSpice\System Identification\FREQ_Cap\Transformer_3Cap_Verification_6outputs.txt", opts);
%dataTable = readtable("C:\Users\q634145\OneDrive - BMW Group\04_Praxiseinsätze\ES-663\04_LTSpice\System Identification\FREQ_Cap\Transformer_3Cap_Verification_6outputs_pwl_150pF.txt", opts);
% Measured Transformer
dataTable = readtable("C:\Users\q634145\OneDrive - BMW Group\04_Praxiseinsätze\ES-663\04_LTSpice\System Identification\MeasuredTransformer\PlanarTransformer_Trans_Final.txt", opts);
t2 = dataTable.time;
Vout_GT = dataTable.Vout;
Vin_GT = dataTable.Vin;
Vw1_GT = dataTable.Vw1;
Iw1_GT = dataTable.Iw1;
Iw2_GT = dataTable.Iw2;
I1_GT = dataTable.I1 * (-1);
I2_GT = dataTable.I2 * (-1);

%%
% Sim and Plot
%V_in_test = 400*sin(2*pi*100000*t2);
%V_in_test = 10*sin(2*pi*50*t2);
V_in_test = Vin_GT;
%V_in_test = 100 * sin(2*pi*100e3 * t2);
% Vout_sim = lsim(ssVo,V_in_test,t2);
% Iw1_sim = lsim(ssIw1,V_in_test,t2);
% Iw2_sim = lsim(ssIw2,V_in_test,t2);
% load vin.mat
% V_in_test = vin.v1;
% t2 = vin.t;

combined_sim = lsim(ssCombined,V_in_test,t2);
Windingcombined_sim = lsim(ssWindingCombined,V_in_test,t2);

figure()
hold on
plot(t2, combined_sim(:,1))
plot(t2, Vout_GT)
legend('Combined Model', 'GT')
title('Vout')

figure()
hold on
plot(t2, combined_sim(:,2))
plot(t2, I2_GT)
legend('Combined Model', 'GT')
title('I1')

figure()
hold on
plot(t2, combined_sim(:,3))
plot(t2, I1_GT)
legend('Combined Model', 'GT')
title('I2')

figure()
hold on
plot(t2, Windingcombined_sim(:,2))
plot(t2, Iw1_GT)
legend('Combined Model', 'GT')
title('Iw1')

figure()
hold on
plot(t2, Windingcombined_sim(:,3))
plot(t2, Iw2_GT)
legend('Combined Model', 'GT')
title('Iw2')

figure()
hold on
plot(t2, Windingcombined_sim(:,1))
plot(t2, Vw1_GT)
legend('Combined Model', 'GT')
title('Vw1')