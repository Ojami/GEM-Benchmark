function varargout = BenchmarkPanel(varargin)
% BenchmarkPanel GUI consists of several benchmarks divided in 2 sections:
% 1- Comparison test, which include gene essentiality, drug response,
% OG/TS/LOF numerator, growth prediction and metabolite uptake/secretion rates
% and 2- Consistency tests, which include resolution power,
% cross-validation, and robustness to noise. Each test is accompanied with
% well-documented functions along with experimental data used to develop
% this panel.
% 
% Oveis Jamialahmad
% TMU. Biotech Dept. 
% Jan 2018
% 
% BENCHMARKPANEL MATLAB code for BenchmarkPanel.fig
%      BENCHMARKPANEL, by itself, creates a new BENCHMARKPANEL or raises the existing
%      singleton*.
%
%      H = BENCHMARKPANEL returns the handle to a new BENCHMARKPANEL or the handle to
%      the existing singleton*.
%
%      BENCHMARKPANEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BENCHMARKPANEL.M with the given input arguments.
%
%      BENCHMARKPANEL('Property','Value',...) creates a new BENCHMARKPANEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BenchmarkPanel_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BenchmarkPanel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BenchmarkPanel

% Last Modified by GUIDE v2.5 18-Sep-2018 21:42:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BenchmarkPanel_OpeningFcn, ...
                   'gui_OutputFcn',  @BenchmarkPanel_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before BenchmarkPanel is made visible.
function BenchmarkPanel_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BenchmarkPanel (see VARARGIN)
clc
currPath = which('BenchmarkPanel.m');
currPath = regexprep(currPath,'BenchmarkPanel.m','');
addpath(genpath([currPath,'General Functions']));

set(handles.figure1,'CloseRequestFcn',@closeGUI)

% Choose default command line output for BenchmarkPanel
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes BenchmarkPanel wait for user response (see UIRESUME)
% uiwait(handles.figure1);
 function closeGUI(hObject, eventdata, handles)
currPath = which('BenchmarkPanel.m');
currPath = regexprep(currPath,'BenchmarkPanel.m','');
rmpath([currPath,'General Functions']);
AT = getappdata(0);
AT = fieldnames(AT);
for i = 1:numel(AT)
    rmappdata(0,AT{i});
end
delete(hObject);


% --- Outputs from this function are returned to the command line.
function varargout = BenchmarkPanel_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in growthPush.
function growthPush_Callback(hObject, eventdata, handles)
% hObject    handle to growthPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
methodChoice = methodSelection(handles);
if ~isempty(methodChoice)
    currPath = which('BenchmarkPanel.m');
    currPath = regexprep(currPath,'BenchmarkPanel.m','');
    addpath(genpath([currPath,'GEMs_Comparison']));
    growthEval(methodChoice,1) % Default model: Recon 1
    rmpath(genpath([currPath,'GEMs_Comparison']));
end



% --- Executes on button press in essentialPush.
function essentialPush_Callback(hObject, eventdata, handles)
% hObject    handle to essentialPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
methodChoice = methodSelection(handles);
if ~isempty(methodChoice)
    currPath = which('BenchmarkPanel.m');
    currPath = regexprep(currPath,'BenchmarkPanel.m','');
    addpath(genpath([currPath,'GEMs_Comparison']));    
    findEssGenes(methodChoice,1) % Default model: Recon 1
    rmpath(genpath([currPath,'GEMs_Comparison']));
end


% --- Executes on button press in exometaPush.
function exometaPush_Callback(hObject, eventdata, handles)
% hObject    handle to exometaPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
methodChoice = methodSelection(handles);
if ~isempty(methodChoice)
    currPath = which('BenchmarkPanel.m');
    currPath = regexprep(currPath,'BenchmarkPanel.m','');
    addpath(genpath([currPath,'GEMs_Comparison']));
    MetChecker(methodChoice,1) % Default model: Recon 1
    rmpath(genpath([currPath,'GEMs_Comparison']));
end


% --- Executes on button press in DrugPush.
function DrugPush_Callback(hObject, eventdata, handles)
% hObject    handle to DrugPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
methodChoice = methodSelection(handles);
if ~isempty(methodChoice)
    currPath = which('BenchmarkPanel.m');
    currPath = regexprep(currPath,'BenchmarkPanel.m','');
    addpath(genpath([currPath,'GEMs_Comparison']));
    datasetList1 = {'HolbeckDataset';'GarnettDataset';'GDSCDataset'};
    datasetList = {'Holbeck et al.';'Garnett et al.';'GDSC'};
    [datasetSel,okButt] = listdlg('PromptString',{'Select target dataset(s):';'For more than one, hold ctrl'},...
        'ListString',datasetList,...
        'Name','Drug Response Dataset','ListSize',[150,100]);
    if ~okButt
        uiwait(errordlg('You should select at least one dataset!','Drug Response Dataset'));
        return
    end
    datasetSel = datasetList1(datasetSel);
    for i = 1:numel(datasetSel)
        fprintf('===================================\n')
        fprintf('Checking for dataset: %s\n',datasetSel{i})
        DrugResponse(datasetSel{i},methodChoice,1) % Default model: Recon 1
        fprintf('===================================\n')
    end
    rmpath(genpath([currPath,'GEMs_Comparison']));
end

% --- Executes on button press in OGPush.
function OGPush_Callback(hObject, eventdata, handles)
% hObject    handle to OGPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
methodChoice = methodSelection(handles);
if ~isempty(methodChoice)
    currPath = which('BenchmarkPanel.m');
    currPath = regexprep(currPath,'BenchmarkPanel.m','');
    addpath(genpath([currPath,'GEMs_Comparison']));
    OncoTsNumerator(methodChoice,1) % Default model: Recon 1
    rmpath(genpath([currPath,'GEMs_Comparison']));
end

% --- Executes on button press in ResPowerPush.
function ResPowerPush_Callback(hObject, eventdata, handles)
% hObject    handle to ResPowerPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
methodChoice = methodSelection(handles);
if ~isempty(methodChoice)
    currPath = which('BenchmarkPanel.m');
    currPath = regexprep(currPath,'BenchmarkPanel.m','');
    addpath(genpath([currPath,'GEMs_Comparison']));
    ResPower(methodChoice)
    rmpath(genpath([currPath,'GEMs_Comparison']));
    fprintf('======= DONE! ========\n')
end


% --- Executes on button press in CVPush.
function CVPush_Callback(hObject, eventdata, handles)
% hObject    handle to CVPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiwait(helpdlg({'To see how to arrange CV models, please';'read "HelponCV.txt" file within GEMs_CV folder'},...
        'CV analysis'));
methodChoice = methodSelection(handles);
if ~isempty(methodChoice)
    currPath = which('BenchmarkPanel.m');
    currPath = regexprep(currPath,'BenchmarkPanel.m','');
    addpath(genpath([currPath,'GEMs_CV']));
    growthEval_CV(methodChoice,1) % Default model: Recon 1
    if any(ismember(methodChoice,'PRIME')) || any(ismember(methodChoice,'TRFBA')) || any(ismember(methodChoice,'pFBA'))
        methodChoice(ismember(methodChoice,'PRIME')) = [];
        methodChoice(ismember(methodChoice,'TRFBA')) = [];
        methodChoice(ismember(methodChoice,'pFBA')) = [];
    end
    CVanalyzer(methodChoice)
    rmpath(genpath([currPath,'GEMs_CV']));
    fprintf('======= DONE! ========\n')
end

% --- Executes on button press in NoisePush.
function NoisePush_Callback(hObject, eventdata, handles)
% hObject    handle to NoisePush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiwait(helpdlg({'To see how to arrange noisy models, please';'read "HelponNoise.txt" file within GEMs_Noisy folder'},...
        'Robustness to noise'));
methodChoice = methodSelection(handles);
if ~isempty(methodChoice)
    currPath = which('BenchmarkPanel.m');
    currPath = regexprep(currPath,'BenchmarkPanel.m','');
    addpath(genpath([currPath,'GEMs_Noisy']));
    growthEval_noisy(methodChoice,1) % Default model: Recon 1
    if any(ismember(methodChoice,'PRIME')) || any(ismember(methodChoice,'TRFBA')) || any(ismember(methodChoice,'pFBA'))
        methodChoice(ismember(methodChoice,'PRIME')) = [];
        methodChoice(ismember(methodChoice,'TRFBA')) = [];
        methodChoice(ismember(methodChoice,'pFBA')) = [];
    end
    ResPower_Noisy(methodChoice)
    rmpath(genpath([currPath,'GEMs_Noisy']));
    fprintf('======= DONE! ========\n')
end

% --- Executes on button press in pFBAcheck.
function pFBAcheck_Callback(hObject, eventdata, handles)
% hObject    handle to pFBAcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pFBAcheck


% --- Executes on button press in TRFBAcheck.
function TRFBAcheck_Callback(hObject, eventdata, handles)
% hObject    handle to TRFBAcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TRFBAcheck


% --- Executes on button press in PRIMEcheck.
function PRIMEcheck_Callback(hObject, eventdata, handles)
% hObject    handle to PRIMEcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PRIMEcheck


% --- Executes on button press in GIMMEcheck.
function GIMMEcheck_Callback(hObject, eventdata, handles)
% hObject    handle to GIMMEcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of GIMMEcheck


% --- Executes on button press in iMATcheck.
function iMATcheck_Callback(hObject, eventdata, handles)
% hObject    handle to iMATcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of iMATcheck


% --- Executes on button press in INITcheck.
function INITcheck_Callback(hObject, eventdata, handles)
% hObject    handle to INITcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of INITcheck


% --- Executes on button press in mCADREcheck.
function mCADREcheck_Callback(hObject, eventdata, handles)
% hObject    handle to mCADREcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mCADREcheck


% --- Executes on button press in FASTCOREcheck.
function FASTCOREcheck_Callback(hObject, eventdata, handles)
% hObject    handle to FASTCOREcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FASTCOREcheck


% --- Executes on button press in FASTCORMICScheck.
function FASTCORMICScheck_Callback(hObject, eventdata, handles)
% hObject    handle to FASTCORMICScheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FASTCORMICScheck


% --- Executes on button press in CORDAcheck.
function CORDAcheck_Callback(hObject, eventdata, handles)
% hObject    handle to CORDAcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CORDAcheck
