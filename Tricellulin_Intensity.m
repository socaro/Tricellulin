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

% Last Modified by GUIDE v2.5 10-Sep-2018 11:36:24

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
%set(handles.axes1,'ActivePositionProperty','position');
handles.pos1=handles.axes1.Position;
handles.pos2=handles.axes2.Position;
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
handles.cellsize=1024/11;
handles.detectjcts=false;
handles.cellrowcurrent=1;
handles.zoom=zoom;
handles.pan=pan;
handles.colors=['g','w','y','m','b','r','c','k','g','w','y','m','b','r','k','g','w','y','m','b','r','c'];
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Tricellulin_Intensity wait for user response (see UIRESUME)
% uiwait(handles.figure1);


function Tricellulin_Intensity_OutputFcn(hObject, eventdata, handles, varargin)
% 
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

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
    case 'q'
        addaux_Callback(hObject,eventdata,handles)
    case 'n'
        newcell_Callback(hObject,eventdata,handles)
    case 'p'
        polygons_Callback(hObject,eventdata,handles)
    case 'x'
        newcellrect_Callback(hObject,eventdata,handles)
    case 'u'
        handles=automatic_cell(handles);
        guidata(hObject,handles);
    case 'c'
        handles=changecellrow(hObject,eventdata,handles);
        guidata(hObject,handles);
    case 'k'
        editpgon_Callback(hObject,eventdata,handles);
    case '`'
        loadsave(hObject,eventdata,handles);
    case 'i'
        handles.cellrowcurrent=handles.cellrowcurrent+1;
        s=handles.cellrowcurrent;
        set(handles.woundrow,'String',num2str(s));
        guidata(hObject,handles);
    case 'm'
        handles=changejunction(handles);
        guidata(hObject,handles);
    case 'o'
        handles.cellrowcurrent=handles.cellrowcurrent-1;
        s=handles.cellrowcurrent;
        set(handles.woundrow,'String',num2str(s));
        guidata(hObject,handles);
    case 'z'
        uitoggletool1_ClickedCallback(hObject, eventdata, handles);
    case 't'
        uitoggletool3_ClickedCallback(hObject, eventdata, handles);
        
        
end

% --- Executes on button press in selectimage.
function selectimage_Callback(hObject, eventdata, handles)
% hObject    handle to selectimage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[imname,pathname]=uigetfile('*.tif','Select B-Catenin Image');
if imname~=0
handles.pathname=pathname;
handles.imname=strrep(imname,'.tif','');
handles.imcat=imread(fullfile(pathname,imname));
handles.imnuc=imread(fullfile(pathname,strrep(imname,'_cat_','_nuc_')));
handles.imtric=imread(fullfile(pathname,strrep(imname,'_cat_','_tric_')));
handles.imfuse=imfuse(imadjust(handles.imcat),imadjust(handles.imnuc),'falsecolor','ColorChannels',[0 1 2]);
%axes(handles.axes1);
xlength=length(handles.imcat(:,1,:));
%pos=get(handles.axes1,'position');
%axis(handles.axes1,[0 pos(3) 0 pos(4)]);
%set(handles.axes1,'position',handles.pos1,'XLim',[0 xlength],'YLim',[0 xlength]);
imshow(handles.imfuse,'Parent',handles.axes1);%,'Parent',handles.axes1,'InitialMagnification',200,'Border','tight','XData',[0 pos(3)],'YData',[0 pos(4)]);%,'InitialMagnification','fit');

%set(handles.axes1,'units','normalized','outerposition',handles.pos1);
%handles.axes1.Position=handles.pos1;
axes(handles.axes2);
%pos=get(handles.axes2,'position');
imshow(imadjust(handles.imtric));
handles.axes2.Position=handles.pos2;
%set(handles.axes2,'position',pos);
handles.cells={};
handles.angles={};
handles.pgons={};
handles.cell_in={};
handles.cell_tric={};
handles.angles_range=[];
handles.cell_tric_avg=[];
handles.cell_njcts=[];
handles.cell_area=[];
handles.cell_perimeter=[];
handles.c=[];
handles.r=[];
handles.truejct=[];
handles.cellrow=[];

%imshowpair(handles.imcat,imadjust(handles.imtric));
set(handles.entersigma,'Enable','on');
set(handles.detecttricellulin,'Enable','on');
set(handles.newcell,'Enable','on');
set(handles.celllist,'Enable','on');
set(handles.imagename,'String',imname);
end
guidata(hObject,handles);

function handles=createcell(handles,in,f)
axes(handles.axes1);
hold on; plot(handles.c(in),handles.r(in),'og');
%handles.pgons{f}=alphaShape(handles.c(in),handles.r(in),'HoleThreshold',1000000);%,'Simplify',true);
set(handles.celllist,'Enable','on');
%in=[in(end); in; in(1)];
x1=handles.c(in);
y1=handles.r(in);
%K=convhull(x,y);
% a=alphaShape(x,y,75);
% [~,P] = boundaryFacets(a);
%handles.pgons{f}=polyshape(P);

k=boundary(x1,y1,0.6);
handles.pgons{f}=polyshape(x1(k),y1(k));
[x,y] = boundary(handles.pgons{f});
if sum(isnan(x))>0
    a=alphaShape(x1,y1,75);
    [~,P] = boundaryFacets(a);
    if length(P)>2
    handles.pgons{f}=polyshape(P);
    [x,y] = boundary(handles.pgons{f});
    if sum(isnan(x))>0
        a=alphaShape(x1,y1,200);
        [~,P] = boundaryFacets(a);
        handles.pgons{f}=polyshape(P);
        [x,y] = boundary(handles.pgons{f});
    end
    else 
        return;
    end
end
[xc,yc]=centroid(handles.pgons{f});c=[xc yc];
x=[x;x(2)];
y=[y;y(2)];
s=handles.s_tric;
handles.cell_in{f}=zeros(1,length(x)-2);
handles.angles{f}=zeros(1,length(x)-2);
handles.cell_tric{f}=zeros(1,length(x)-2);
truejct_c=false(1,length(x)-2);
for i=1:length(x)-2
   c1=[x(i) y(i)];c2=[x(i+2) y(i+2)];c3=[x(i+1) y(i+1)];
   r31=c1-c3;r32=c2-c3;r3c=c-c3;
   %handles.angles{f}(i)=acos(dot(r31,r32)/(norm(r31)*norm(r32))); 
   a13c=acos(dot(r31,r3c)/(norm(r31)*norm(r3c)));
   ac32=acos(dot(r32,r3c)/(norm(r32)*norm(r3c)));
   a=a13c+ac32;
   handles.angles{f}(i)=a;
   handles.cell_tric{f}(i)=mean(mean(handles.imtric(y(i+1)-s:y(i+1)+s,x(i+1)-s:x(i+1)+s)));
   innew=knnsearch([handles.c(in) handles.r(in)],[x(i+1) y(i+1)]);
   innew=in(innew);
   if handles.truejct(innew)
       truejct_c(i)=true;
   end
   handles.cell_in{f}(i)=innew;
end
%handles.cell_in{f}=handles.cell_in{f}(truejct_c);
handles.angles{f}=handles.angles{f}(truejct_c);
handles.cell_tric{f}=handles.cell_tric{f}(truejct_c);
handles.cell_njcts(f)=length(handles.cell_in{f}(truejct_c));
handles.cell_tric_avg(f)=mean(handles.cell_tric{f});
handles.angles_range(f)=range(handles.angles{f})*(180/pi);
handles.cell_area(f)=area(handles.pgons{f});
handles.cell_perimeter(f)=perimeter(handles.pgons{f});
handles.cellrow(f)=handles.cellrowcurrent;
hold on; plot(handles.pgons{f},'FaceColor',handles.colors(handles.cellrow(f)),'FaceAlpha',0.3,'EdgeAlpha',0.5);
axes(handles.axes2);
hold on; plot(handles.pgons{f},'FaceColor',handles.colors(handles.cellrow(f)),'FaceAlpha',0.3,'EdgeAlpha',0.5);
set(handles.info,'String',sprintf('Cell: %d, Tricellulin Intensity: %d, Avg all cells: %d', f,handles.cell_tric_avg(f),mean(handles.cell_tric_avg)));
set(handles.polygons,'Enable','on');
set(handles.clearjcts,'Enable','off');



function handles = changecellrow(hObject,eventdata,handles)
axes(handles.axes1)
f=0;
[c,r]=my_ginput;
for n=1:length(c)
for i=1:length(handles.cell_in)
    if isinterior(handles.pgons{i},c(n),r(n))
        f=i;
        break;
    end
end
handles.cells{f}=sprintf('Cell %d',f');
set(handles.celllist,'String',handles.cells);
set(handles.celllist,'Value',f);
handles.cellrow(f)=handles.cellrowcurrent;
end
polygons_Callback(hObject,eventdata,handles);

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


function handles=automatic_cell(handles)
    
    axes(handles.axes1);
    h=imrect;
    
    set(handles.deletejunctions,'Enable','off');
    
    [in,jct,handles.truejct]=celledges_tric(handles.imcat,[handles.c handles.r],handles.cellsize,false,h,handles.detectjcts,handles.truejct,handles.imtric);
    handles.c=jct(:,1);
    handles.r=jct(:,2);
    
    axes(handles.axes1);
    cla(handles.axes1);
    imshow(handles.imfuse);
    hold on; scatter(handles.c,handles.r,200,'.','r')
    axes(handles.axes2);
    cla(handles.axes2);
    imshow(imadjust(handles.imtric));
    hold on; scatter(handles.c,handles.r,100,'.','b');
    
    
    for i=1:length(in)
        if length(in{i})>2
            f=length(handles.cells)+1;
            handles.cells{f}=sprintf('Cell %d',f');
            set(handles.celllist,'String',handles.cells);
            set(handles.celllist,'Value',f);
            handles=createcell(handles,in{i},f);
        end
    end
    delete(h);


% --- Executes on button press in newcell.
function newcellrect_Callback(hObject, eventdata, handles)
% hObject    handle to newcell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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
jctangles=cell(length(handles.c),1);
cellrow=handles.cellrow;

in=cell(size(cell_jcts));
angles_range=handles.angles_range;
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
        polyout = scale(handles.pgons{i},1.1,[cx cy]);
        [x,y]=boundary(polyout);
        BW=poly2mask(x,y,size(handles.imtric,1),size(handles.imtric,2));
        tric_mask=uint16(BW).*handles.imtric.*uint16(bwcat);
        tric_avg=mean(mean(tric_mask(tric_mask~=0)));
        signal{i}=cell_tric{i}(1:length(angles{i}));
        signal{i}=signal{i}./tric_avg;
        cell_tric_avg_n(i)=handles.cell_tric_avg(i)/tric_avg;
        in{i}=cell_jcts{i}(handles.truejct(cell_jcts{i})==1);
        for j=1:length(in{i})
            jctangles{in{i}(j)}=[jctangles{in{i}(j)} angles{i}(j)];
        end
    end
    cell_tric_avg_abs=handles.cell_tric_avg;
end
anglestd=zeros(1,length(jctangles));
anglerange=zeros(1,length(jctangles));
for i=1:length(jctangles)
    if abs(sum(jctangles{i})-2*pi)<0.1
    anglestd(i)=std(jctangles{i}.*(180/pi));
    anglerange(i)=range(jctangles{i}.*(180/pi));
    end
    
end
%signalangle=cell(length(handles.r),1);
anglesmat=zeros(length(cells),10);
% signalmat=zeros(length(cells),10);
anglestd_jct=cell(1,length(cells));
anglestd_range=cell(1,length(cells));
signal_std=zeros(length(cells),1);
for i=1:length(cells)
   if length(angles{i})>10
       anglesmat(i,1:10)=[angles{i}(1:10)].*(180/pi);
   else
   anglesmat(i,1:length(angles{i}))=[angles{i}].*(180/pi);
   end
%    signalmat(i,1:length(cell_tric{i}))=[cell_tric{i}];
   anglestd_jct{i}=anglestd(in{i});
   anglerange_jct{i}=anglerange(in{i});
   signal_std(i)=std(signal{i});
end
anglescol=cell2mat(angles).*(180/pi);
signalcol=cell2mat(signal);

anglestdcol=cell2mat(anglestd_jct);
anglerangecol=cell2mat(anglerange_jct);
umtopix=handles.umtopix;
cell_area=handles.cell_area*umtopix^2;
cell_perimeter=handles.cell_perimeter*umtopix;
cell_sf=cell_perimeter./(cell_area.^0.5);
x=handles.c;
y=handles.r;
cell_njcts=handles.cell_njcts;
pgons=handles.pgons;
truejct=handles.truejct;
save(fullfile(handles.pathname,name),'cells','cell_area','cell_perimeter','cell_njcts','pgons','cell_jcts','cell_tric','cell_tric_avg_n','cell_tric_avg_abs','x','y','umtopix','imname','angles','angles_range','truejct','anglescol','anglestdcol','anglerangecol','signalcol','signal_std','cellrow');
A=[{name,'','','','','','','','','','','','','','','','','','';'cell_njcts','cell_tric_avg','cell_tric_avg_n','signal_std','cell_area','cell_perimeter','shape_factor','angles_range','cellrow','angle1','angle2','angle3','angle4','angle5','angle6','angle7','angle8','angle9','angle10'};mat2cell([cell_njcts(1:length(cells))',cell_tric_avg_abs(1:length(cells))',cell_tric_avg_n(1:length(cells))',signal_std(1:length(cells)),cell_area(1:length(cells))',cell_perimeter(1:length(cells))',cell_sf(1:length(cells))',angles_range(1:length(cells))',cellrow(1:length(cells))',anglesmat],ones(length(cells),1),ones(1,19))];
pathname=which('Tricellulin_Intensity.m');
pathname=pathname(1:end-23);
excelname=[pathname 'Excels/TricellulinIntensity' imname(1:8) '.xlsx'];
type=imname(end-15:end-5);
% if exist(excelname,'file')==2
%     %excel=xlsread('TricellulinIntensity.xlsx');
%     [~,sheets,~] = xlsfinfo(excelname);
%     %A={excel;name,'';'cell_njcts','cell_tric_avg';cell_njcts,cell_tric_avg,};
%     xlswrite(excelname,A,length(sheets)+1);
% else
    xlswrite(excelname,A,type);
% end

ind=find(anglestdcol);
%Asignal=[{name,'','','','','','','','','','','','','','','','','','','';'angle1','angle2','angle3','angle4','angle5','angle6','angle7','angle8','angle9','angle10','signal1','signal2','signal3','signal4','signal5','signal6','signal7','signal8','signal9','signal10'};mat2cell([anglesmat,signalmat],ones(length(cell_njcts),1),ones(1,20))];
Asignal=[{name,'','','';'angle','angle std','angle range','signal'};mat2cell([anglescol(ind)',anglestdcol(ind)',anglerangecol(ind)',signalcol(ind)'],ones(length(signalcol(ind)),1),ones(1,4))];

excelname=[pathname 'Excels\TricellulinAngle' imname(1:8) '.xlsx'];

% if exist(excelname,'file')==2
%     %excel=xlsread('TricellulinIntensity.xlsx');
%     [~,sheets,~] = xlsfinfo(excelname);
%     %A={excel;name,'';'cell_njcts','cell_tric_avg';cell_njcts,cell_tric_avg,};
%     xlswrite(excelname,A,length(sheets)+1);
% else
    xlswrite(excelname,Asignal,type);
% end

%xlswrite('SignalAngle.xlsx',signalangle);

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
% imbw=imbinarize(handles.imcat(:,:,1),'adaptive','Sensitivity',0.4);
% imbw=bwareaopen(imbw,100,4);
% [ycoords, xcoords, ~, ~] = findendsjunctions(imbw,0);

handles.c=xcoords;
handles.r=ycoords;
handles.truejct=true(size(handles.r));
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

function handles=selectpgon(handles)
axes(handles.axes1)
[c,r]=ginput(1);
for i=1:length(handles.cell_in)
    if inpolygon(c,r,handles.c(handles.cell_in{i}),handles.r(handles.cell_in{i}))==1
        f=i;
        handles.cells{f}=sprintf('Cell %d',f');
        set(handles.celllist,'String',handles.cells);
        set(handles.celllist,'Value',f);
        h=impoly;
        Pos=getPosition(h);
        in=find(inpolygon(handles.c,handles.r,Pos(:,1),Pos(:,2)));
        handles=createcell(handles,in,f);
        break;
    end
end




% --- Executes on button press in addjunctions.
function addjunctions_Callback(hObject, eventdata, handles)
% hObject    handle to addjunctions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
[c,r]=my_ginput;

	cn=[];
    rn=[];
    s=5;
    for i=1:length(c)
        [cc, rr]=nodesmax(handles.imtric(r(i)-s:r(i)+s,c(i)-s:c(i)+s),handles.addsigma,2);
        if isempty(cc)
            cn=[cn; c(i)];
            rn=[rn; r(i)];
        else
            cn=[cn; cc'+c(i)-s];
            rn=[rn; rr'+r(i)-s];
        end
    end
%[c,r,~]=impixel(handles.im);
handles.c=[handles.c; cn];
handles.r=[handles.r; rn];
handles.truejct=[handles.truejct; true(size(cn))];
%imshow(handles.im);
 axes(handles.axes1);
 cla(handles.axes1);
 imshow(handles.imfuse);
 hold on; scatter(handles.c(handles.truejct==1),handles.r(handles.truejct==1),200,'.','r');
 hold on; scatter(handles.c(~handles.truejct==1),handles.r(~handles.truejct==1),200,'.','b');
 axes(handles.axes2);
 cla(handles.axes2);
 imshow(imadjust(handles.imtric));
 %imshowpair(handles.imcat,imadjust(handles.imtric));
 hold on; scatter(handles.c(handles.truejct==1),handles.r(handles.truejct==1),100,'.','b');
 hold on; scatter(handles.c(~handles.truejct==1),handles.r(~handles.truejct==1),100,'.','g');
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
    [IDX,~] = knnsearch([handles.c,handles.r],[x,y]);
    delete_nodes=zeros(size(handles.c));
    delete_nodes(IDX)=1;
    handles.c=handles.c(~delete_nodes);
    handles.r=handles.r(~delete_nodes);
    handles.truejct=handles.truejct(~delete_nodes);
    axes(handles.axes1);
    cla(handles.axes1);
    imshow(handles.imfuse);
    hold on; scatter(handles.c(handles.truejct==1),handles.r(handles.truejct==1),200,'.','r');
    hold on; scatter(handles.c(~handles.truejct==1),handles.r(~handles.truejct==1),200,'.','b');
    axes(handles.axes2);
    cla(handles.axes2);
    imshow(imadjust(handles.imtric));
    %imshowpair(handles.imcat,imadjust(handles.imtric));
     hold on; scatter(handles.c(handles.truejct==1),handles.r(handles.truejct==1),100,'.','b');
 hold on; scatter(handles.c(~handles.truejct==1),handles.r(~handles.truejct==1),100,'.','g');
 end
 guidata(hObject,handles);

    function handles = changejunction(handles)
       [x,y] = ginput(1);%('Color','r');
       [xnew,ynew]=ginput(1);
       if ~isempty(x)&&~isempty(xnew)
        [IDX,~] = knnsearch([handles.c,handles.r],[x,y]);
        handles.c(IDX)=xnew;handles.r(IDX)=ynew;
        for i=1:length(handles.cell_in)
            if ~isempty(find(handles.cell_in{i}==IDX))
            handles=createcell(handles,handles.cell_in{i},i);
            end
        end
       end

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
set(gcf,'Pointer','arrow');
axes(handles.axes1);
cla(handles.axes1);
imshow(handles.imfuse);
hold on; scatter(handles.c(handles.truejct==1),handles.r(handles.truejct==1),200,'.','r');
 hold on; scatter(handles.c(~handles.truejct==1),handles.r(~handles.truejct==1),200,'.','b');
for i=1:length(handles.cells)
    hold on; plot(handles.pgons{i},'FaceColor',handles.colors(handles.cellrow(i)),'FaceAlpha',0.3,'EdgeAlpha',0.5);
end
axes(handles.axes2);
cla(handles.axes2);
imshow(imadjust(handles.imtric));
 hold on; scatter(handles.c(handles.truejct==1),handles.r(handles.truejct==1),100,'.','b');
 hold on; scatter(handles.c(~handles.truejct==1),handles.r(~handles.truejct==1),100,'.','g');
for i=1:length(handles.cells)
    hold on; plot(handles.pgons{i},'FaceColor',handles.colors(handles.cellrow(i)),'FaceAlpha',0.3,'EdgeAlpha',0.5);
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


function handles=loaddata(handles,matname,pathname)
if matname~=0
load(fullfile(pathname,matname),'imname','cells','cell_area','cell_tric_avg_abs','cell_perimeter','cell_njcts','pgons','cell_jcts','cell_tric','x','y','umtopix','angles','angles_range');
try
  load(fullfile(pathname,matname),'truejct');
  handles.truejct=truejct;
catch
  handles.truejct=true(size(x));
end
try
    load(fullfile(pathname,matname),'cellrow');
    handles.cellrow=cellrow;
catch
    handles.cellrow=ones(size(cell_njcts));
end
handles.imname=strrep(imname,'cat_*','cat_');
imname=[imname '.tif'];
set(handles.imagename,'String',imname);
handles.cell_tric_avg=cell_tric_avg_abs;
handles.angles=angles;
handles.cellrowcurrent=1;
%handles.angles_range=angles_range;
handles.cells=cells;
handles.pgons=pgons;
handles.cell_njcts=cell_njcts;
handles.umtopix=umtopix;
for i=1:length(cells)
handles.cell_area(i)=area(handles.pgons{i});
handles.cell_perimeter(i)=perimeter(handles.pgons{i});
handles.angles_range(i)=range(handles.angles{i})*(180/pi);
end
%handles.cell_area=cell_area/umtopix;
%handles.cell_perimeter=cell_perimeter/umtopix;
handles.pgons=pgons;
handles.c=x;
handles.r=y;
handles.cell_in=cell_jcts;
handles.cell_tric=cell_tric;
handles.pathname=pathname;
handles.imcat=imread(fullfile(pathname,imname));
handles.imnuc=imread(fullfile(pathname,strrep(imname,'_cat_','_nuc_')));
handles.imtric=imread(fullfile(pathname,strrep(imname,'_cat_','_tric_')));
handles.imfuse=imfuse(imadjust(handles.imcat),imadjust(handles.imnuc),'falsecolor','ColorChannels',[0 1 2]);

axes(handles.axes1);
imshow(handles.imfuse);
hold on; scatter(handles.c(handles.truejct==1),handles.r(handles.truejct==1),200,'.','r');
 hold on; scatter(handles.c(~handles.truejct==1),handles.r(~handles.truejct==1),200,'.','b');
axes(handles.axes2);
imshow(imadjust(handles.imtric));
 hold on; scatter(handles.c(handles.truejct==1),handles.r(handles.truejct==1),100,'.','b');
 hold on; scatter(handles.c(~handles.truejct==1),handles.r(~handles.truejct==1),100,'.','g');
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
end


% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[matname,pathname]=uigetfile('*.mat','Select recently saved data');
handles=loaddata(handles,matname,pathname);
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
handles.cell_tric={};
handles.angles_range=[];
handles.cell_tric_avg=[];
handles.cell_njcts=[];
handles.cell_area=[];
handles.cell_perimeter=[];

axes(handles.axes1);
imshow(handles.imfuse);
hold on; scatter(handles.c(handles.truejct==1),handles.r(handles.truejct==1),200,'.','r');
 hold on; scatter(handles.c(~handles.truejct==1),handles.r(~handles.truejct==1),200,'.','b');
 
axes(handles.axes2);
imshow(imadjust(handles.imtric));
 hold on; scatter(handles.c(handles.truejct==1),handles.r(handles.truejct==1),100,'.','b');
 hold on; scatter(handles.c(~handles.truejct==1),handles.r(~handles.truejct==1),100,'.','g');
set(handles.celllist,'String',handles.cells);
set(handles.addjunctions,'Enable','on');
set(handles.deletejunctions,'Enable','on');
set(handles.polygons,'Enable','off');
set(handles.clearjcts,'Enable','on');
guidata(hObject,handles);







% --- Executes on button press in addaux.
function addaux_Callback(hObject, eventdata, handles)
% hObject    handle to addaux (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
hold on;
[cn,rn]=my_ginput;

%[c,r,~]=impixel(handles.im);
handles.c=[handles.c; cn];
handles.r=[handles.r; rn];
handles.truejct=[handles.truejct; false(size(cn))];
%imshow(handles.im);
 axes(handles.axes1);
 cla(handles.axes1);
 imshow(handles.imfuse);
 hold on; scatter(handles.c(handles.truejct==1),handles.r(handles.truejct==1),200,'.','r');
 hold on; scatter(handles.c(~handles.truejct==1),handles.r(~handles.truejct==1),200,'.','b');
 axes(handles.axes2);
 cla(handles.axes2);
 imshow(imadjust(handles.imtric));
 %imshowpair(handles.imcat,imadjust(handles.imtric));
 hold on; scatter(handles.c(handles.truejct==1),handles.r(handles.truejct==1),100,'.','b');
 hold on; scatter(handles.c(~handles.truejct==1),handles.r(~handles.truejct==1),100,'.','g');
guidata(hObject,handles);


% --- Executes on button press in editpgon.
function editpgon_Callback(hObject, eventdata, handles)
% hObject    handle to editpgon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1)
f=0;
while f<1
[c,r]=ginput(1);
for i=1:length(handles.pgons)
    if isinterior(handles.pgons{i},c,r)
        f=i;
        break;
    end
end
end
handles.cells{f}=sprintf('Cell %d',f');
set(handles.celllist,'String',handles.cells);
set(handles.celllist,'Value',f);
h=impoly;
Pos=getPosition(h);
in=find(inpolygon(handles.c,handles.r,Pos(:,1),Pos(:,2)));
handles=createcell(handles,in,f);
delete(h);
guidata(hObject,handles);



function loadsave(hObject,eventdata,handles)
[matname,pathname]=uigetfile('*.mat','Select recently saved data','MultiSelect','on');
for i=1:length(matname)
    handles=loaddata(handles,matname{i},pathname);
    save_Callback(hObject,eventdata,handles);
end


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes on button press in clearjcts.
function clearjcts_Callback(hObject, eventdata, handles)
% hObject    handle to clearjcts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.c=[];
handles.r=[];
handles.truejct=[];
axes(handles.axes1);
    cla(handles.axes1);
    imshow(handles.imfuse);
    hold on; scatter(handles.c,handles.r,200,'.','r')
    axes(handles.axes2);
    cla(handles.axes2);
    imshow(imadjust(handles.imtric));
    hold on; scatter(handles.c,handles.r,100,'.','b');
    guidata(hObject,handles);
    



function cellnr_Callback(hObject, eventdata, handles)
% hObject    handle to cellnr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
s=str2double(get(hObject,'String'));
if ~isnan(s)
handles.cellsize=round(length(handles.imcat(:,:,1))/s);
end
guidata(hObject,handles);
% Hints: get(hObject,'String') returns contents of cellnr as text
%        str2double(get(hObject,'String')) returns contents of cellnr as a double


% --- Executes during object creation, after setting all properties.
function cellnr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cellnr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in autojct.
function autojct_Callback(hObject, eventdata, handles)
% hObject    handle to autojct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.detectjcts=get(hObject,'Value');
guidata(hObject,handles);
% Hint: get(hObject,'Value') returns toggle state of autojct



function woundrow_Callback(hObject, eventdata, handles)
% hObject    handle to woundrow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
s=str2double(get(hObject,'String'));
if ~isnan(s)
handles.cellrowcurrent=s;
end
guidata(hObject,handles);
% Hints: get(hObject,'String') returns contents of woundrow as text
%        str2double(get(hObject,'String')) returns contents of woundrow as a double


% --- Executes during object creation, after setting all properties.
function woundrow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to woundrow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --------------------------------------------------------------------
function uitoggletool1_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% if strcmp(handles.zoom.Enable,'on')
%     handles.zoom.Enable='off';
% else
    handles.zoom.Enable='on';
    handles.uitoggletool1.Visible='off';
    handles.uitoggletool3.Visible='off';
    while true
    k=waitforbuttonpress;
    if k==1
    handles.zoom.Enable='off';
    handles.uitoggletool1.Visible='on';
    handles.uitoggletool3.Visible='on';
    break;
    end
    end
% end


% --------------------------------------------------------------------
function uitoggletool3_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% if strcmp(handles.pan.Enable,'on')
%     handles.pan.Enable='off';
% else

    handles.pan.Enable='on';
    handles.uitoggletool3.Visible='off';
    handles.uitoggletool1.Visible='off';
    while true
    k=waitforbuttonpress;
    if k==1
    handles.pan.Enable='off';
    handles.uitoggletool3.Visible='on';
    handles.uitoggletool1.Visible='on';
    break;
    end
    end
% end
