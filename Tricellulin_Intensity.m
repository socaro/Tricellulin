function varargout = Tricellulin_Intensity(varargin)
%TRICELLULIN_INTENSITY MATLAB code file for Tricellulin_Intensity.fig
%      TRICELLULIN_INTENSITY, by itself, creates a new TRICELLULIN_INTENSITY or raises the existing
%      singleton*.
%
%      H = TRICELLULIN_INTENSITY returns the handle to a new TRICELLULIN_INTENSITY or the handle to
%      the existing singleton*.
%
%      TRICELLULIN_INTENSITY('Property','Value',...) creates a new TRICELLULIN_INTENSITY using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to Tricellulin_Intensity_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      TRICELLULIN_INTENSITY('CALLBACK') and TRICELLULIN_INTENSITY('CALLBACK',hObject,...) call the
%      local function named CALLBACK in TRICELLULIN_INTENSITY.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Tricellulin_Intensity

% Last Modified by GUIDE v2.5 02-May-2018 14:53:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Tricellulin_Intensity_OpeningFcn, ...
                   'gui_OutputFcn',  @Tricellulin_Intensity_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before Tricellulin_Intensity is made visible.
function Tricellulin_Intensity_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for Tricellulin_Intensity

handles.output = hObject;
set(handles.entersigma,'Enable','off');
set(handles.detecttricellulin,'Enable','off');
set(handles.newcell,'Enable','off');
set(handles.addjunctions,'Enable','off');
set(handles.deletejunctions,'Enable','off');
set(handles.save,'Enable','off');
set(handles.celllist,'Enable','off');
set(handles.polygons,'Enable','off');
set(handles.objective,'String',[{'40x Objective'};{'60x Objective'}]);
handles.umtopix=0.325;
handles.s_tric=1;
handles.addsigma=1.5;
handles.binarytric=false;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Tricellulin_Intensity wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Tricellulin_Intensity_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
switch eventdata.Character
    case 'd'
        deletejunctions_Callback(hObject, eventdata, handles)
    case 'a'
        addjunctions_Callback(hObject, eventdata, handles)
    case 'n'
        newcell_Callback(hObject,eventdata,handles)
    case 'p'
        polygons_Callback(hObject,eventdata,handles)
    case 'x'
        newcellrect_Callback(hObject,eventdata,handles)
end

% --- Executes on button press in selectimage.
function selectimage_Callback(hObject, eventdata, handles)
% hObject    handle to selectimage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[imname,pathname]=uigetfile('*.tif','Select B-Catenin Image');
handles.pathname=pathname;
handles.imname=strrep(imname,'.tif','');
handles.imcat=imread(fullfile(pathname,imname));
handles.imnuc=imread(fullfile(pathname,strrep(imname,'_cat_','_nuc_')));
handles.imtric=imread(fullfile(pathname,strrep(imname,'_cat_','_tric_')));
handles.imfuse=imfuse(imadjust(handles.imcat),imadjust(handles.imnuc),'falsecolor','ColorChannels',[0 1 2]);
axes(handles.axes1);
imshow(handles.imfuse);
axes(handles.axes2);
imshow(imadjust(handles.imtric));
%imshowpair(handles.imcat,imadjust(handles.imtric));
set(handles.entersigma,'Enable','on');
set(handles.detecttricellulin,'Enable','on');
set(handles.newcell,'Enable','on');
set(handles.celllist,'Enable','on');
guidata(hObject,handles);

function handles=createcell(handles,in,f)
axes(handles.axes1);
hold on; plot(handles.c(in),handles.r(in),'og');
handles.cell_in{f}=in;
%handles.pgons{f}=alphaShape(handles.c(in),handles.r(in),'HoleThreshold',1000000);%,'Simplify',true);
set(handles.celllist,'Enable','on');
s=handles.s_tric;
for i=1:length(in)
    handles.cell_tric{f}(i)=mean(mean(handles.imtric(handles.r(in(i))-s:handles.r(in(i))+s,handles.c(in(i))-s:handles.c(in(i))+s)));
end
%in=[in(end); in; in(1)];


handles.cell_tric_avg(f)=mean(handles.cell_tric{f});
handles.cell_njcts(f)=length(handles.cell_in{f});
x=handles.c(in);
y=handles.r(in);
%K=convhull(x,y);
a=alphaShape(x,y,50);
[~,P] = boundaryFacets(a);
handles.pgons{f}=polyshape(P);
[x,y] = boundary(handles.pgons{f});
[xc,yc]=centroid(handles.pgons{f});c=[xc yc];
x=[x;x(2)];
y=[y;y(2)];
for i=1:length(x)-2
   c1=[x(i) y(i)];c2=[x(i+2) y(i+2)];c3=[x(i+1) y(i+1)];
   r31=c1-c3;r32=c2-c3;r3c=c-c3;
   %handles.angles{f}(i)=acos(dot(r31,r32)/(norm(r31)*norm(r32))); 
   a13c=acos(dot(r31,r3c)/(norm(r31)*norm(r3c)));
   ac32=acos(dot(r32,r3c)/(norm(r32)*norm(r3c)));
   handles.angles{f}(i)=a13c+ac32;
end
handles.angles_range(f)=range(handles.angles{f});
handles.cell_area(f)=area(handles.pgons{f});
handles.cell_perimeter(f)=perimeter(handles.pgons{f});
hold on; plot(handles.pgons{f},'FaceColor','r','FaceAlpha',0.3,'EdgeAlpha',0.5);
set(handles.info,'String',sprintf('Cell: %d, Tricellulin Intensity: %d, Avg all cells: %d', f,handles.cell_tric_avg(f),mean(handles.cell_tric_avg)));
set(handles.polygons,'Enable','on');


% --- Executes on button press in newcell.
function newcell_Callback(hObject, eventdata, handles)
% hObject    handle to newcell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.deletejunctions,'Enable','off');
f=length(handles.cells)+1;
handles.cells{f}=sprintf('Cell %d',f');
set(handles.celllist,'String',handles.cells);
set(handles.celllist,'Value',f);
% axes(handles.axes1);
% imshow(handles.im);
% hold on; plot(handles.c,handles.r,'og');
% for i=1:length(handles.cells)-1
%     hold on; plot(handles.pgons{i},'FaceColor','g','FaceAlpha',0.3,'EdgeAlpha',0.5);
% end
axes(handles.axes1);
h=impoly;
Pos=getPosition(h);
in=find(inpolygon(handles.c,handles.r,Pos(:,1),Pos(:,2)));
handles=createcell(handles,in,f);
delete(h);
guidata(hObject,handles);

% --- Executes on button press in newcell.
function newcellrect_Callback(hObject, eventdata, handles)
% hObject    handle to newcell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.deletejunctions,'Enable','off');
f=length(handles.cells)+1;
handles.cells{f}=sprintf('Cell %d',f');
set(handles.celllist,'String',handles.cells);
set(handles.celllist,'Value',f);
% axes(handles.axes1);
% imshow(handles.im);
% hold on; plot(handles.c,handles.r,'og');
% for i=1:length(handles.cells)-1
%     hold on; plot(handles.pgons{i},'FaceColor','g','FaceAlpha',0.3,'EdgeAlpha',0.5);
% end
axes(handles.axes1);
h=imrect;
Pos=getPosition(h);
in=find(inpolygon(handles.c,handles.r,[Pos(1) Pos(1) Pos(1)+Pos(3) Pos(1)+Pos(3) Pos(1)],[Pos(2) Pos(2)+Pos(4) Pos(2)+Pos(4) Pos(2) Pos(2)]));
handles=createcell(handles,in,f);
delete(h);
guidata(hObject,handles);

% --- Executes on selection change in celllist.
function celllist_Callback(hObject, eventdata, handles)
% hObject    handle to celllist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f=get(hObject,'Value');
axes(handles.axes1);
% imshowpair(handles.imcat,imadjust(handles.imtric));
% hold on; plot(handles.c(handles.cell_in{f}),handles.r(handles.cell_in{f}),'.r');
% for i=1:length(handles.cells)
%     hold on; plot(handles.pgons{i},'FaceColor','g','FaceAlpha',0.3,'EdgeAlpha',0.5);
% end
%hold on; plot(handles.pgons{f},'FaceColor','r','FaceAlpha',0.3,'EdgeAlpha',0.5);
h=impoly;
Pos=getPosition(h);
in=find(inpolygon(handles.c,handles.r,Pos(:,1),Pos(:,2)));
handles=createcell(handles,in,f);
delete(h);
guidata(hObject,handles);
% Hints: contents = cellstr(get(hObject,'String')) returns celllist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from celllist


% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imname=handles.imname;
cells=handles.cells;
cell_jcts=handles.cell_in;
cell_tric=handles.cell_tric;
angles=handles.angles;
anglesmat=zeros(length(cells),10);
for i=1:length(cells)
   anglesmat(i,1:length(angles{i}))=[angles{i}].*(180/pi);
end
angles_range=handles.angles_range.*(180/pi);
bwcat=imbinarize(handles.imcat,'adaptive');
%tric_avg=mean(mean(handles.imtric));
if handles.binarytric==true
    name=[handles.imname handles.nameaddon 'binary.mat'];
    for i=1:length(cell_jcts)
        %BW=roipoly(handles.imtric,handles.c(cell_jcts{i}),handles.r(cell_jcts{i}));
%         [x,y]=boundary(handles.pgons{i});
%         BW=poly2mask(x,y,size(handles.imtric,1),size(handles.imtric,2));
%         tric_mask=uint16(BW).*handles.imtric.*uint16(bwcat);
%         tric_avg=mean(mean(tric_mask(tric_mask~=0)));
        [cx,cy]=centroid(handles.pgons{i});
        polyout = scale(handles.pgons{i},1.4,[cx cy]);
        [x,y]=boundary(polyout);
        BW=poly2mask(x,y,size(handles.imtric,1),size(handles.imtric,2));
        tric_mask=uint16(BW).*handles.imtric.*uint16(bwcat);
        tric_avg=mean(mean(tric_mask(tric_mask~=0)));
        cell_tric_avg_n(i)=mean(handles.cell_tric{i}>1.3*tric_avg);
    end
    cell_tric_avg_abs=handles.cell_tric_avg;
else
    name=[handles.imname handles.nameaddon '.mat'];
    for i=1:length(cell_jcts)
        [cx,cy]=centroid(handles.pgons{i});
        polyout = scale(handles.pgons{i},1.4,[cx cy]);
        [x,y]=boundary(polyout);
        BW=poly2mask(x,y,size(handles.imtric,1),size(handles.imtric,2));
        tric_mask=uint16(BW).*handles.imtric.*uint16(bwcat);
        tric_avg=mean(mean(tric_mask(tric_mask~=0)));
        cell_tric_avg_n(i)=handles.cell_tric_avg(i)/tric_avg;
    end
    cell_tric_avg_abs=handles.cell_tric_avg;
end
umtopix=handles.umtopix;
cell_area=handles.cell_area*umtopix^2;
cell_perimeter=handles.cell_perimeter*umtopix;
cell_sf=cell_perimeter./(cell_area.^0.5);
x=handles.c;
y=handles.r;
cell_njcts=handles.cell_njcts;
pgons=handles.pgons;
save(fullfile(handles.pathname,name),'cells','cell_area','cell_perimeter','cell_njcts','pgons','cell_jcts','cell_tric','cell_tric_avg_n','cell_tric_avg_abs','x','y','umtopix','imname','angles','angles_range');
if exist('TricellulinIntensity.xlsx','file')==2
    %excel=xlsread('TricellulinIntensity.xlsx');
    [~,sheets,~] = xlsfinfo('TricellulinIntensity.xlsx');
    %A={excel;name,'';'cell_njcts','cell_tric_avg';cell_njcts,cell_tric_avg,};
    A=[{name,'','','','','','','','','','','','','','','','';'cell_njcts','cell_tric_avg','cell_tric_avg_n','cell_area','cell_perimeter','shape_factor','angles_range','angle1','angle2','angle3','angle4','angle5','angle6','angle7','angle8','angle9','angle10'};mat2cell([cell_njcts',cell_tric_avg_abs',cell_tric_avg_n',cell_area',cell_perimeter',cell_sf',angles_range',anglesmat],ones(length(cell_njcts),1),ones(1,17))];
    xlswrite('TricellulinIntensity.xlsx',A,length(sheets)+1);
else
    A=[{name,'','','','','','','','','','','','','','','','';'cell_njcts','cell_tric_avg','cell_tric_avg_n','cell_area','cell_perimeter','shape_factor','angles_range','angle1','angle2','angle3','angle4','angle5','angle6','angle7','angle8','angle9','angle10'};mat2cell([cell_njcts',cell_tric_avg_abs',cell_tric_avg_n',cell_area',cell_perimeter',cell_sf',angles_range',anglesmat],ones(length(cell_njcts),1),ones(1,17))];
    xlswrite('TricellulinIntensity.xlsx',A);
end

guidata(hObject,handles);


function entersigma_Callback(hObject, eventdata, handles)
% hObject    handle to entersigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.sigma=str2double(get(hObject,'String'));
set(handles.detecttricellulin,'Enable','on');
guidata(hObject,handles);
% Hints: get(hObject,'String') returns contents of entersigma as text
%        str2double(get(hObject,'String')) returns contents of entersigma as a double


% --- Executes during object creation, after setting all properties.
function entersigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to entersigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in detecttricellulin.
function detecttricellulin_Callback(hObject, eventdata, handles)
% hObject    handle to detecttricellulin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.cells=cell(0);
[xcoords,ycoords]=nodesmax(handles.imtric,handles.sigma,3);
handles.c=xcoords;
handles.r=ycoords;
handles.cell_in=cell(0);

axes(handles.axes1);
cla(handles.axes1);
imshow(handles.imfuse);
hold on; scatter(xcoords,ycoords,200,'.','r')
axes(handles.axes2);
cla(handles.axes2);
imshow(imadjust(handles.imtric));
%imshowpair(handles.imcat,imadjust(handles.imtric));
hold on; scatter(xcoords,ycoords,100,'.','b');

set(handles.deletejunctions,'Enable','on');
guidata(hObject,handles);

% --- Executes on button press in addjunctions.
function addjunctions_Callback(hObject, eventdata, handles)
% hObject    handle to addjunctions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
hold on;
[c,r]=my_ginput;
if handles.addsigma==0
    cn=c;
    rn=r;
else
	cn=[];
    rn=[];
    s=10;
    for i=1:length(c)
        [cc, rr]=nodes(handles.imtric(r(i)-s:r(i)+s,c(i)-s:c(i)+s),handles.addsigma);
        cn=[cn; cc+c(i)-s];
        rn=[rn; rr+r(i)-s];
    end
end
%[c,r,~]=impixel(handles.im);
handles.c=[handles.c; cn];
handles.r=[handles.r; rn];
%imshow(handles.im);
 axes(handles.axes1);
 cla(handles.axes1);
 imshow(handles.imfuse);
 hold on; scatter(handles.c,handles.r,200,'.','r');
 axes(handles.axes2);
 cla(handles.axes2);
 imshow(imadjust(handles.imtric));
 %imshowpair(handles.imcat,imadjust(handles.imtric));
 hold on; scatter(handles.c,handles.r,100,'.','b');
guidata(hObject,handles);


% --- Executes on button press in deletejunctions.
function deletejunctions_Callback(hObject, eventdata, handles)
% hObject    handle to deletejunctions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
 %zoom off;zoom off;pan off;datacursormode off;
 [x,y] = my_ginput;%('Color','r');
 if ~isempty(x)
 [IDX,D] = knnsearch([handles.c,handles.r],[x,y]);
 delete_nodes=zeros(size(handles.c));
 delete_nodes(IDX)=1;
 handles.c=handles.c(~delete_nodes);
 handles.r=handles.r(~delete_nodes);
 axes(handles.axes1);
 cla(handles.axes1);
 imshow(handles.imfuse);
 hold on; scatter(handles.c,handles.r,200,'.','r');
 axes(handles.axes2);
 cla(handles.axes2);
 imshow(imadjust(handles.imtric));
 %imshowpair(handles.imcat,imadjust(handles.imtric));
 hold on; scatter(handles.c,handles.r,100,'.','b');
 end
 guidata(hObject,handles);



function addjctsigma_Callback(hObject, eventdata, handles)
% hObject    handle to addjctsigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.addsigma=str2double(get(hObject,'String'));
set(handles.addjunctions,'Enable','on');
guidata(hObject,handles);
% Hints: get(hObject,'String') returns contents of addjctsigma as text
%        str2double(get(hObject,'String')) returns contents of addjctsigma as a double


% --- Executes during object creation, after setting all properties.
function addjctsigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to addjctsigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function name_Callback(hObject, eventdata, handles)
% hObject    handle to name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.nameaddon=get(hObject,'String');
set(handles.save,'Enable','on');
guidata(hObject,handles);
% Hints: get(hObject,'String') returns contents of name as text
%        str2double(get(hObject,'String')) returns contents of name as a double


% --- Executes during object creation, after setting all properties.
function name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in polygons.
function polygons_Callback(hObject, eventdata, handles)
% hObject    handle to polygons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
cla(handles.axes1);
imshow(handles.imfuse);
hold on; scatter(handles.c,handles.r,200,'.','r');
for i=1:length(handles.cells)
    hold on; plot(handles.pgons{i},'FaceColor','g','FaceAlpha',0.3,'EdgeAlpha',0.5);
end


% --- Executes on button press in binary.
function binary_Callback(hObject, eventdata, handles)
% hObject    handle to binary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.binarytric=get(hObject,'Value');
guidata(hObject,handles);

% if handles.binarytric==false
%     handles.binarytric=true;
% else
%     handles.binarytric=false;
% end

% Hint: get(hObject,'Value') returns toggle state of binary


% --- Executes on selection change in objective.
function objective_Callback(hObject, eventdata, handles)
% hObject    handle to objective (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
obj=get(hObject,'Value');
switch obj
    case 1
        handles.umtopix=0.325;
    case 2
        handles.umtopix=0.2167;     
end
guidata(hObject,handles);
        
% Hints: contents = cellstr(get(hObject,'String')) returns objective contents as cell array
%        contents{get(hObject,'Value')} returns selected item from objective


% --- Executes during object creation, after setting all properties.
function objective_CreateFcn(hObject, eventdata, handles)
% hObject    handle to objective (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[matname,pathname]=uigetfile('*.mat','Select recently saved data');
load(fullfile(pathname,matname),'imname','cells','cell_area','cell_perimeter','cell_njcts','pgons','cell_jcts','cell_tric','x','y','umtopix','angles','angles_range');
imname=[imname '.tif'];
handles.angles=angles;
handles.angles_range=angles_range;
handles.cells=cells;
handles.pgons=pgons;
handles.cell_njcts=cell_njcts;
handles.umtopix=umtopix;
handles.cell_area=cell_area/umtopix;
handles.cell_perimeter=cell_perimeter/umtopix;
handles.pgons=pgons;
handles.c=x;
handles.r=y;
handles.cell_in=cell_jcts;
handles.cell_tric=cell_tric;
handles.pathname=pathname;
handles.imname=imname;
handles.imcat=imread(fullfile(pathname,imname));
handles.imnuc=imread(fullfile(pathname,strrep(imname,'_cat_','_nuc_')));
handles.imtric=imread(fullfile(pathname,strrep(imname,'_cat_','_tric_')));
handles.imfuse=imfuse(imadjust(handles.imcat),imadjust(handles.imnuc),'falsecolor','ColorChannels',[0 1 2]);

axes(handles.axes1);
imshow(handles.imfuse);
hold on; scatter(handles.c,handles.r,200,'.','r');
axes(handles.axes2);
imshow(imadjust(handles.imtric));
hold on; scatter(handles.c,handles.r,100,'.','b');
%imshowpair(handles.imcat,imadjust(handles.imtric));
if ~isempty(handles.cells)
set(handles.deletejunctions,'Enable','off');
set(handles.polygons,'Enable','on');
f=length(handles.cells);
set(handles.celllist,'String',handles.cells);
set(handles.celllist,'Value',f);
else
set(handles.entersigma,'Enable','on');
set(handles.detecttricellulin,'Enable','on');
end
set(handles.newcell,'Enable','on');
set(handles.celllist,'Enable','on');
guidata(hObject,handles);


% --- Executes on button press in clear.
function clear_Callback(hObject, eventdata, handles)
% hObject    handle to clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.cells={};
handles.angles={};
handles.pgons={};
handles.cell_in={};
handles.angles_range=[];
handles.cell_tric_avg=[];
handles.cell_njcts=[];
handles.cell_area=[];
handles.cell_perimeter=[];

axes(handles.axes1);
imshow(handles.imfuse);
hold on; scatter(handles.c,handles.r,200,'.','r');
set(handles.celllist,'String',handles.cells);
set(handles.addjunctions,'Enable','on');
set(handles.deletejunctions,'Enable','on');
set(handles.polygons,'Enable','off');
guidata(hObject,handles);
