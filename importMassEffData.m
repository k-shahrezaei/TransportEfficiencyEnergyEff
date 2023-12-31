%% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: /Users/k.shahrezaei/Documents/MATLAB/TransportEfficiencyEnergyEff/massEffData.xlsx
%    Worksheet: Blad1
%
% Auto-generated by MATLAB on 18-Aug-2023 15:39:15

%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 9);

% Specify sheet and range
opts.Sheet = "Blad1";
opts.DataRange = "A4:I14";

% Specify column names and types
opts.VariableNames = ["TruckType", "gwMin", "gwMax", "ewMin", "ewMax", "payloadCapMax", "payloadCapRatio", "massEffMIn", "massEfmax"];
opts.VariableTypes = ["categorical", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, "TruckType", "EmptyFieldRule", "auto");

% Import the data
massEffData = readtable("/Users/k.shahrezaei/Documents/MATLAB/TransportEfficiencyEnergyEff/massEffData.xlsx", opts, "UseExcel", false);


%% Clear temporary variables
clear opts