function varargout = ARNquant(varargin)
% ARNQUANT MATLAB code for ARNquant.fig
%      ARNQUANT, by itself, creates a new ARNQUANT or raises the existing
%      singleton*.
%
%      H = ARNQUANT returns the handle to a new ARNQUANT or the handle to
%      the existing singleton*.
%
%      ARNQUANT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ARNQUANT.M with the given input arguments.
%
%      ARNQUANT('Property','Value',...) creates a new ARNQUANT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ARNquant_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ARNquant_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ARNquant

% Last Modified by GUIDE v2.5 08-Feb-2016 19:18:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ARNquant_OpeningFcn, ...
                   'gui_OutputFcn',  @ARNquant_OutputFcn, ...
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

end
% --- Executes just before ARNquant is made visible.
function ARNquant_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ARNquant (see VARARGIN)

% Choose default command line output for ARNquant
handles.output = hObject;
zoom(handles.spotax, 'off')



% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ARNquant wait for user response (see UIRESUME)
% uiwait(handles.figure1);

end
% --- Outputs from this function are returned to the command line.
function varargout = ARNquant_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

% --- Executes on button press in apply.
function apply_Callback(hObject, eventdata, handles)
% hObject    handle to apply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles, 'rna') && ~isempty(handles.rna)

    rna = setrnanuc(handles.rna, handles.maskfile);
    result = [];
    switch get(get(handles.algochoice,'SelectedObject'),'Tag')
          case 'mean',  result = mean_algo(handles, rna);
          case 'median',  result = median_algo(handles, rna);
          case 'classification',  result = class_algo(handles, rna);
    end
    if ~isempty(result)
        handles.result = double(result);
        guidata(hObject, handles);
    end
end

end

function rna  = class_algo(handles, rna)
background = 0;
ind = get(handles.algorithm,'Value');
nclass = get(handles.centroid, 'Value');
if nclass < 3
    nclass = nclass + 1;
else
    nclass = 0;
end

intensity = rna(:,3);
cytospot = get(handles.cytospot, 'Value');
nucspot = get(handles.nucspot, 'Value');
coeff = str2double(get(handles.coeff, 'String'));
if isnan(coeff)
    coeff=1.5;
end

if cytospot && ~nucspot
    intensity = rna(rna(:,4)==background, 3);
elseif ~cytospot && nucspot
    intensity = rna(rna(:,4)>background, 3);
end

test_class = [2, 3, 4, 5];
sil_score = zeros(4, 1);

if nclass == 0
    for i = 1:numel(test_class)
        try 
            cluster = get_cluster(intensity, ind, test_class(i));
            sil_score(i) = mean(silhouette(intensity, cluster));
        catch me
            disp(me)
        end
    end
    [max_sil, sil_pos] = max(sil_score);
    nclass = test_class(sil_pos);
end

cluster = get_cluster(intensity, ind, nclass);
cluster_uniq = sort(unique(cluster));
ncluster =  numel(cluster_uniq);
cluster_uniq(:,2) = zeros(size(cluster_uniq));

hf = figure;
hp = uipanel('Units', 'normalized', 'Position',[0 0.2 1 0.8], 'Parent', hf);
subplot(2,1,1, 'Parent', hp);
cc = summer(ncluster);
clustleg = {};
for i = 1:ncluster
    cluster_uniq(i, 2) = mean(intensity(cluster == cluster_uniq(i,1))); % set mean intensity for each cluster
    clust_int_dist = intensity(cluster == cluster_uniq(i,1));
    binsize = str2double(get(handles.binsize, 'String'));
    if isnan(binsize)
        binsize = 500;
    end
    binranges = min(clust_int_dist):binsize:max(clust_int_dist)+100;
    [bincounts] = histc(clust_int_dist,binranges);
    bar(binranges,bincounts, 'FaceColor', cc(i, :), 'EdgeColor', 'none');
    clustleg{i} = strcat(num2str(i), ' cluster');
    hold on;
end

xlabel('Culster distribution','FontSize',12,'FontWeight','bold');
ylabel('size','FontSize',12,'FontWeight','bold');
legend(clustleg);

[~, guessed_single] = min(abs(cluster_uniq(:,2) - median(intensity)));
mean_int = cluster_uniq(guessed_single, 2);

rna(:,end+1)=(double(rna(:,3)/mean_int)>coeff) & (rna(:,4)~=background);
rna(:,end+1)=round(rna(:,3)/mean_int);

subplot(2,1,2,'Parent', hp);
nascents = sort(unique(rna(:,end)))';
cc = hsv(numel(nascents));
leg = {};

binsize = str2double(get(handles.binsize, 'String'));
if isnan(binsize)
    binsize = 500;
end
for j = 1:numel(nascents)
    intensity = rna(rna(:,end) == nascents(j), 3);
    binranges = min(intensity):binsize:max(intensity)+100;
    [bincounts] = histc(intensity,binranges);
    bh = bar(binranges,bincounts, 'FaceColor', cc(j, :), 'EdgeColor', 'none');
    leg{j} = strcat(num2str(nascents(j)), ' nascents');
    hold on;
end

yl = ylim;
plt = plot([mean_int,mean_int], yl, 'k-', 'LineWidth',2);
xlabel('Intensity','FontSize',12,'FontWeight','bold');
ylabel('Count','FontSize',12,'FontWeight','bold');
legend(leg);
title('Intensity distribution between spots')

control = uipanel('Units', 'normalized','Position', [0 0 1 0.15], 'Background', 'white');

bg = uibuttongroup('Title','Best cluster for single intensity', 'Parent', control);
for i = 1:length(cluster_uniq(:,1))
    rb = uicontrol('Style','radiobutton', 'Parent', bg, 'String', num2str(i), 'Units', 'normalized', 'Position', [0.1*i-0.05 0.5 0.1 0.35]);
    if i == guessed_single
        set(rb, 'Value', (get(rb, 'Max'))); 
    end
end

set(bg,'SelectionChangeFcn',{@updateVal});

function updateVal(~, ~)
    choice = str2double(get(get(bg,'SelectedObject'),'String'));
    mean_int = cluster_uniq(choice,2);
    set(plt,'XData', [mean_int, mean_int]);
end

hbut=uicontrol('Style', 'pushbutton', 'String', 'Ok', 'Parent',control, 'Units', 'normalized', 'Position', [0.85 0.08 0.1 0.4], 'Callback', {@accept, hf});

function accept(~, ~, fig)
    close(fig);
end

waitfor(hbut, 'UserData')
rna(:,end-1)=(double(rna(:,3)/mean_int)>coeff) & (rna(:,4)~=background);
rna(:,end)=round(rna(:,3)/mean_int);

end

function cluster = get_cluster(intensity, ind, nclass)
cluster = zeros(size(intensity));
switch ind
    case 1, cluster = ckmeans_algo(intensity, nclass); % ckmeans
    case 2, cluster = gmm_algo(intensity, nclass);% gaussian mixture
    case 3, cluster = kmeans_algo(intensity, nclass); % kmeans
    case 4, cluster = kde_algo(intensity, nclass); % kde

end
end

function cluster = ckmeans_algo(intensity, nclass)
cluster = ckmeans(intensity, nclass);

end

function cluster = gmm_algo(intensity, nclass)
gmodel = gmdistribution.fit(intensity, nclass);
cluster = gmodel.cluster(intensity);
end

function cluster = kmeans_algo(intensity, nclass)
cluster = kmeans(intensity, nclass);

end

function rna = median_algo(handles, rna)

background = 0;
deviate = str2double(get(handles.std, 'String'));
if isnan(deviate)
    deviate = 2;
end
outliers = abs(rna(:,3) - median(rna(:,3))) > deviate*std(rna(:,3));
cytospot = get(handles.cytospot, 'Value');
nucspot = get(handles.nucspot, 'Value');

s_intent = median(rna(~outliers,3));
if cytospot && ~nucspot
    s_intent = median(rna(rna(:,4)==background & ~outliers, 3));
elseif ~cytospot && nucspot
    s_intent = median(rna(rna(:,4)>background & ~outliers, 3));
    % if you don't check any of the box
    % which is strange, it will just ignore your choice
end
    
coeff = str2double(get(handles.coeff, 'String'));
if isnan(coeff)
    coeff=1.5;
end
rna = intdistmean(rna, s_intent, coeff);

end

function rna = mean_algo(handles, rna)
min_inten =  get(handles.minslider, 'Value');
max_inten = get(handles.maxslider, 'Value');
background = 0; % confident on using that value as the background
formean_rna = rna(rna(:,3)>=min_inten & rna(:, 3)<=max_inten, :);
s_intent = str2double(get(handles.meanint, 'String'));
if isnan(s_intent)
    cytospot = get(handles.cytospot, 'Value');
    nucspot = get(handles.nucspot, 'Value');
    if cytospot && ~nucspot
        formean_rna = rna(rna(:,4)==background,  :);
    elseif ~cytospot && nucspot
        formean_rna = rna(rna(:,4)>background, :);
        % if you don't check any of the box
        % which is strange, it will just ignore your choice
    end
    
    s_intent = mean(formean_rna(:,3));
end
coeff = str2double(get(handles.coeff, 'String'));
if isnan(coeff)
    coeff=1.5;
end
rna = intdistmean(rna, s_intent, coeff);

end

function rna = intdistmean(rna, mean_int, coeff)
% only nuclear spot can be set as trasncription site
% rna : X Y intent nucleus istrans nascentnumber
background = 0;
rna(:,end+1)=(double(rna(:,3)/mean_int)>coeff) & (rna(:,4)~=background);
rna(:,end+1)=round(rna(:,3)/mean_int);

hf = figure;
%hist(intensity);
nascents = sort(unique(rna(:,end)))';
cc = hsv(numel(nascents));
leg = {};
binsize = str2double(get(handles.binsize, 'String'));
if isnan(binsize)
    binsize = 500;
end
for j=1:numel(nascents)
    intensity = rna(rna(:,end) == nascents(j), 3);
    binranges = min(intensity):binsize:max(intensity)+100;
    [bincounts] = histc(intensity,binranges);
    bh = bar(binranges,bincounts, 'FaceColor', cc(j, :), 'EdgeColor', 'none');
    leg{j} = strcat(num2str(nascents(j)), ' nascents');
    hold on;
end
yl = ylim;
disp(mean_int)
plot([mean_int,mean_int], yl, 'k-', 'LineWidth',2);
xlabel('Intensity','FontSize',12,'FontWeight','bold');
ylabel('Count','FontSize',12,'FontWeight','bold');
legend(leg);
title('Intensity distribution between spots')
end

function coor=setrnanuc(coor, mask)
% Trouver le noyau de chaque spot
nuc_int = int64(sort(unique(mask(mask>0))));
nuc_vec = zeros(numel(coor(:,1)),1);
for i=1:length(coor(:,1))
    i_nuc=mask(round(coor(i,2)),round(coor(i,1)));
    proche_nuc=zeros(2,2);
    if(coor(i,2)>1 && coor(i,1)>1)
        proche_nuc= mask(round(coor(i,2))-1:round(coor(i,2))+1,round(coor(i,1))-1:round(coor(i,1))+1);
    end
    if(i_nuc>0)
        nuc_vec(i) = find(nuc_int==i_nuc);
    elseif sum(proche_nuc(:))>0
        nuc_vec(i) = find(nuc_int==proche_nuc(find(proche_nuc~=0,1)));
    end
end
coor(:, end+1)= nuc_vec';
end

% --- Executes on button press in ok.
function ok_Callback(hObject, eventdata, handles)
% hObject    handle to ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles, 'result')
    [filename, pathname] = uiputfile('intensity.locx', 'Save Intensity file');
    dlmwrite(fullfile(pathname, filename), handles.result,'delimiter', '\t', '-append');
    if isfield(handles, 'maskfile') && ~isempty(handles.maskfile)
       saveimage(handles.maskfile, fullfile(pathname, 'label.eps')) 
    end
else
   wdlg = warndlg('Nothing to save, please press Apply before', 'Empty result'); 
   uiwait(wdlg)
end
end

% --- Executes on button press in cancel.
function cancel_Callback(hObject, eventdata, handles)
% hObject    handle to cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg('Are you sure you want to close this program?',...
    'Close Request','Yes','No','Yes');
switch selection
    case 'Yes',
        close(gcf);
    case 'No'
        return
end %switch

end

function std_Callback(hObject, eventdata, handles)
% hObject    handle to std (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of std as text
%        str2double(get(hObject,'String')) returns contents of std as a double

end

% --- Executes during object creation, after setting all properties.
function std_CreateFcn(hObject, eventdata, handles)
% hObject    handle to std (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

% --- Executes on selection change in algorithm.
function algorithm_Callback(hObject, eventdata, handles)
% hObject    handle to algorithm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns algorithm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from algorithm
end

% --- Executes during object creation, after setting all properties.
function algorithm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to algorithm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on selection change in centroid.
function centroid_Callback(hObject, eventdata, handles)
% hObject    handle to centroid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns centroid contents as cell array
%        contents{get(hObject,'Value')} returns selected item from centroid

end

% --- Executes during object creation, after setting all properties.
function centroid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to centroid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

% --- Executes on slider movement.
function minslider_Callback(hObject, eventdata, handles)
% hObject    handle to minslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
currVal = double(get(hObject, 'Value'));
set(handles.minval, 'String', num2str(currVal));

if isfield(handles, 'minplot') 
    lineplot(handles.minplot, currVal)
end
reduce_spot(handles, currVal, get(handles.maxslider, 'Value'))
end

% --- Executes during object creation, after setting all properties.
function minslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

end

% --- Executes on slider movement.
function maxslider_Callback(hObject, eventdata, handles)
% hObject    handle to maxslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
currVal = double(get(hObject, 'Value'));
set(handles.maxval, 'String', num2str(currVal));
if isfield(handles, 'maxplot') 
    lineplot(handles.maxplot, currVal)
end
reduce_spot(handles, get(handles.minslider, 'Value'), currVal);
end

% --- Executes during object creation, after setting all properties.
function maxslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end

function minval_Callback(hObject, eventdata, handles)
% hObject    handle to minval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minval as text
%        str2double(get(hObject,'String')) returns contents of minval as a double

% only accept value in range min to max
currVal = str2double(get(hObject, 'String'));

max_val = get(handles.maxslider, 'Value');
min_val = get(handles.minslider, 'Min');

% oh god, why didn't I use if clauses ??
currVal = (currVal >= min_val && currVal < max_val)*currVal  + (currVal >= max_val)*max_val + (currVal < min_val)*min_val;
set(hObject, 'String', num2str(currVal));
set(handles.minslider, 'Value', currVal)
if isfield(handles, 'minplot') 
    lineplot(handles.minplot, currVal)
end
reduce_spot(handles, currVal, max_val);

end

% --- Executes during object creation, after setting all properties.
function minval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function maxval_Callback(hObject, eventdata, handles)
% hObject    handle to maxval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxval as text
%        str2double(get(hObject,'String')) returns contents of maxval as a double

currVal = str2double(get(hObject, 'String'));
max_val = get(handles.maxslider, 'Max');
min_val = get(handles.minslider, 'Value');
currVal = (currVal > min_val && currVal <= max_val)*currVal  + (currVal > max_val)*max_val + (currVal <= min_val)*min_val;
set(hObject, 'String', num2str(currVal));
set(handles.maxslider, 'Value', currVal)
if isfield(handles, 'maxplot') 
    lineplot(handles.maxplot, currVal)
end
reduce_spot(handles, min_val, currVal);

end

% --- Executes during object creation, after setting all properties.
function maxval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end


function meanint_Callback(hObject, eventdata, handles)
% hObject    handle to meanint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of meanint as text
%        str2double(get(hObject,'String')) returns contents of meanint as a double
end

% --- Executes during object creation, after setting all properties.
function meanint_CreateFcn(hObject, eventdata, handles)
% hObject    handle to meanint (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

% --------------------------------------------------------------------
function mask_Callback(hObject, eventdata, handles)
% hObject    handle to mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname]= uigetfile({'*.tif;*.tiff','Microscopy image file'; '*.*','All (*.*)'}, 'Pick the mask file');

if ~ (isequal(filename,0) || isequal(pathname,0))
    mask = fullfile(pathname, filename);
    if  ishandle(handles.maskax) 
        set(handles.filename, 'String', mask);
        finfo = imfinfo(mask);
        mask = imread(mask);
        if ~ (strcmp(finfo.ColorType, 'grayscale') && min(unique(mask))==0)
            mask = bwlabel(im2bw(mat2gray(mask),0),4);
        end
        handles.maskfile = mask;
        axes(handles.maskax);
        imshow(handles.maskfile, [], 'InitialMagnification', 'fit')
    end
    guidata(hObject, handles);
end

end

function replot(handles, binsize)
    [min_int, min_ind] = min(handles.rna(:,3));
    [max_int, max_ind] = max(handles.rna(:,3));
    axes(handles.spotax);
    binranges = min_int:binsize:max_int+100;
    [bincounts] = histc(handles.rna(:,3),binranges);
    bar(binranges,bincounts,'w');
end

% --------------------------------------------------------------------
function spot_Callback(hObject, eventdata, handles)
% hObject    handle to spot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[rnafile, pathname]= uigetfile({'*.loc;*locx','Localize file (.loc)'; '*.*','All (*.*)'}, 'Pick the mRNA file');
% data should be in the following format:  x y intensity
if ~ (isequal(rnafile,0) || isequal(pathname,0))
    rnaf = fullfile(pathname, rnafile);
    set(handles.spottext, 'String', rnaf);
    rna = load(rnaf);
    cla(handles.spotax);
    handles.rna = rna(:,1:3);
    [min_int, min_ind] =min(rna(:,3));
    [max_int, max_ind] =max(rna(:,3));
    axes(handles.spotax);
    binsize = str2double(get(handles.binsize, 'String'));
    if isnan(binsize)
        binsize = 500;
    end
    binranges = min_int:binsize:max_int+100;
    [bincounts] = histc(rna(:,3),binranges);
    bar(binranges,bincounts,'w');
    
    yl = ylim;
    set(handles.minslider, 'Min', min_int , 'Max', max_int, 'Value', min_int, 'SliderStep',[100 1000]./(max_int-min_int));
    set(handles.maxslider, 'Min', min_int , 'Max', max_int, 'Value', max_int, 'SliderStep',[100 1000]./(max_int-min_int));
    set(handles.minval, 'String', num2str(min_int));
    set(handles.maxval, 'String', num2str(max_int));
    hold on;
    handles.minplot = plot([min_int,min_int], yl, 'g-', 'LineWidth',1.5);
    hold on;
    handles.maxplot = plot([max_int,max_int], yl, 'r-', 'LineWidth',1.5);
    if isfield(handles, 'maskfile') && numel(handles.maskfile) > 0
        handles.spotmask = zeros(size(handles.maskfile));
        axes(handles.maskax)
        hold on;
        handles.spotplot = plot(int64(handles.rna(:,1)),int64(handles.rna(:,2)),'*','Color',[1 0.5 0.5],'MarkerSize',3);
        axes(handles.spotax)
    
    else
        max_size = int64(max(max(rna(:, 1:2))));
        handles.spotmask = zeros(max_size, max_size);
        
    end
    linear_pos = sub2ind(size(handles.spotmask), int32(rna(:,2)), int32(rna(:,1)));
    handles.spotmask(linear_pos) = rna(:,3);   

    guidata(hObject, handles);
end

end


% --------------------------------------------------------------------
function new_Callback(hObject, eventdata, handles)
% hObject    handle to new (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ARNquant

end

function lineplot(h, xval)
if ishandle(h)
    set(h,'XData',[xval xval]);
end

end

function reduce_spot(handles, min_inten, max_inten)
 if isfield(handles, 'spotplot') && ishandle(handles.maskax) 
     rnalist = handles.rna(handles.rna(:,3)>=min_inten & handles.rna(:,3)<=max_inten, :);
     set(handles.spotplot,'XData',int64(rnalist(:,1)))
     set(handles.spotplot,'YData', int64(rnalist(:,2)))
 end
 
end


 % --------------------------------------------------------------------
function uitoggletool2_OnCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h = zoom;
setAllowAxesZoom(h, handles.spotax, 0);
guidata(hObject, handles)
end

% --------------------------------------------------------------------
function uitoggletool3_OnCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = zoom;
setAllowAxesZoom(h, handles.spotax, 0);
guidata(hObject, handles)
end

% --------------------------------------------------------------------
function uitoggletool1_OnCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.dcm_obj = datacursormode(gcf);
set(handles.dcm_obj,'Enable','on', 'UpdateFcn',{@myupdatefcn,hObject});
guidata(hObject, handles)
end

function txt = myupdatefcn(~, event_obj, hFigure)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).
txt = '';
hAxesParent  = get(get(event_obj,'Target'),'Parent');
% get the handles structure for the figure/GUI
handles = guidata(hFigure);
% which axes are we in?
pos = get(event_obj,'Position');

if hAxesParent == handles.maskax    
    txt = {['X: ',num2str(pos(1))],...
           ['Y: ',num2str(pos(2))]};
    intnt = handles.spotmask(pos(2), pos(1));
    if isfield(handles, 'spotmask') && intnt
         txt = {['X: ',num2str(pos(1))],...
           ['Y: ',num2str(pos(2))], ...
          ['Intensity: ',num2str(intnt)]};
    end
else
    % do nothing 
end
end

% --------------------------------------------------------------------
function uitoggletool1_OffCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.dcm_obj,'Enable','off')

% We are going to remove all object that are hggroup but not bar plot.
c = findall(gcf,'Type','hggroup');
bars = findobj('-property','Barwidth');
delete(setdiff(c, bars))
set(gcf,'Pointer','arrow');
end

function saveimage(nucmask, outputfile)
cell_list= sort(unique(nucmask(:)))';
cell_list = cell_list(cell_list~=0);
centroid_list = zeros(length(cell_list),2);

for i=1:length(cell_list)
    bw = nucmask == cell_list(i);
    stat = regionprops(bw, 'Centroid');
    centroid_list(i,:)= stat(1).Centroid;
end

f = figure('color','white','units','normalized','position',[.1 .1 .8 .8]);
bwim =  nucmask <= 0;
imshow(bwim);
set(f,'units','pixels','position',[0 0 size(bwim,1)  size(bwim,2)],'visible','off')
axis off

for i=1:length(cell_list)
    text('position',centroid_list(i,:),'FontWeight','bold','fontsize',10,'string',num2str(i), 'color', [0.5,0.5,0.5]) ;
end
print(f,'-depsc','-r150',outputfile);
close(f);
end

% --- Executes during object creation, after setting all properties.
function coeff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to coeff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in nucspot.
function nucspot_Callback(hObject, eventdata, handles)
% hObject    handle to nucspot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nucspot
end

% --- Executes on button press in cytospot.
function cytospot_Callback(hObject, eventdata, handles)
% hObject    handle to cytospot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cytospot
end


% --- Executes during object creation, after setting all properties.
function spotax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spotax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate spotax
end


% --- Executes during object creation, after setting all properties.
function maskax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maskax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate maskax
end



function binsize_Callback(hObject, eventdata, handles)
% hObject    handle to binsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of binsize as text
%        str2double(get(hObject,'String')) returns contents of binsize as a double

binsize = str2double(get(hObject,'String'));
if isnan(binsize)
    binsize = 500;
end
replot(handles, binsize);
guidata(hObject, handles)

end
% --- Executes during object creation, after setting all properties.
function binsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to binsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
