function varargout = NADH(varargin)
% NADHaxes MATLAB code for NADHaxes.fig
%      NADHaxes, by itself, creates a new NADHaxes or raises the existing
%      singleton*.
%
%      H = NADHaxes returns the handle to a new NADHaxes or the handle to
%      the existing singleton*.
%
%      NADHaxes('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NADHaxes.M with the given input arguments.
%
%      NADHaxes('Property','Value',...) creates a new NADHaxes or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NADH_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NADH_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NADHaxes

% Last Modified by GUIDE v2.5 16-Aug-2019 13:12:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @NADH_OpeningFcn, ...
    'gui_OutputFcn',  @NADH_OutputFcn, ...
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


% --- Executes just before NADHaxes is made visible.
function NADH_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NADHaxes (see VARARGIN)
global black CROPIMAGE_check DETECTNecrosis_check ...
    DETECTViableTHRESH_check SAVE_check FolderSelect_check DETECTViable_check ...
    DetectNecrosisTHRESH_check black2 CreateStruct FontSIZE DATA
DATA={};
black=zeros(3000,3000);
axes(handles.OriginalImageCroppedImage); %set to show on axes1
imshow(black); %display image
axes(handles.NADHaxes); %set to show on axes1
imshow(black); %display image
axes(handles.NecroticFraction); %set to show on axes1
imshow(black); %display image
axes(handles.ViableBorder); %set to show on axes1
imshow(black); %display image
axes(handles.FolderName); %set to show on axes1
insert=ones(5,40);
black2=cat(3, insert*0,insert*0,insert*.7);
imshow(black2); %display image
axes(handles.NecroticFractionAxes);
imshow(black2);
axes(handles.NecroticFraction);
imshow(black)
axes(handles.NecroticBorder);
imshow(black)
CROPIMAGE_check=0;
DETECTViable_check=0;
DETECTNecrosis_check=0;
SAVE_check=0;
FolderSelect_check=0;
DETECTViableTHRESH_check=0;
DetectNecrosisTHRESH_check=0;
FontSIZE=15;
CreateStruct.Interpreter = 'tex';
CreateStruct.WindowStyle = 'modal';
% Choose default command line output for NADHaxes
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes NADHaxes wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = NADH_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in CropImage.
function CropImage_Callback(hObject, eventdata, handles)
global original image_crop black position CROPIMAGE_check FolderSelect_check ...
    CreateStruct FontSIZE
if FolderSelect_check==0;
    msgbox(['\fontsize{' ,num2str(FontSIZE),'}Please select a folder.'], 'Warning','warn',CreateStruct);
elseif FolderSelect_check==1;
    CROPIMAGE_check=1;
    % hObject    handle to CropImage (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    axes(handles.OriginalImageCroppedImage); %set to show on axes1
    h=imrect(gca, [450 1 500 500]); %create a rectangle to get coordinates
    position=wait(h); %get the position of the rectangle
    image_crop=imcrop(original,position); %crop the image
    axes(handles.OriginalImageCroppedImage); %set to show on axes1
    imshow(image_crop); %display image
    axes(handles.NADHaxes); %set to show on axes1
    imshow(black); %display image
    axes(handles.NecroticFraction); %set to show on axes1
    imshow(black); %display image
    axes(handles.ViableBorder); %set to show on axes1
    imshow(black); %display image
end





% --- Executes on button press in Align.
function Align_Callback(hObject, eventdata, handles)
% hObject    handle to Align (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global original CROPIMAGE_check  image_crop FolderSelect_check ...
    CreateStruct FontSIZE
if FolderSelect_check==0;
    msgbox(['\fontsize{' ,num2str(FontSIZE),'}Please select a folder.'], 'Warning','warn',CreateStruct);
elseif FolderSelect_check==1;
    if CROPIMAGE_check==1;%call global image
        [xI,yI]=ginput(2); %gather click location information
        Hdistimage=abs(xI(2)-xI(1)); %calculate horizontal distance
        Vdistimage=abs(yI(1)-yI(2)); %calculate vertical distance
        if Hdistimage>Vdistimage %if horizontal is longer than vertical
            opposite=yI(2)-yI(1); %calculate length of the opposite side
            adjacent=xI(2)-xI(1); %calculate length of the adjacent side
            angle=atand(opposite/adjacent); %calcualte the angle
            image_crop=imrotate(image_crop,angle,'bilinear','crop'); %rotate the image
        else %otherwise
            adjacent=yI(2)-yI(1); %calculate lenght of the opposite side
            opposite=xI(2)-xI(1); %calculate length of the adjacent side
            angle=atand(opposite/adjacent); %calcualte the angle
            image_crop=imrotate(image_crop,-angle,'bilinear','crop'); %rotate the image
        end
        axes(handles.OriginalImageCroppedImage); %set to show on axes1
        imshow(image_crop); %display image
    else
        [xI,yI]=ginput(2); %gather click location information
        Hdistimage=abs(xI(2)-xI(1)); %calculate horizontal distance
        Vdistimage=abs(yI(1)-yI(2)); %calculate vertical distance
        if Hdistimage>Vdistimage %if horizontal is longer than vertical
            opposite=yI(2)-yI(1); %calculate length of the opposite side
            adjacent=xI(2)-xI(1); %calculate length of the adjacent side
            angle=atand(opposite/adjacent); %calcualte the angle
            original=imrotate(original,angle,'bilinear','crop'); %rotate the image
        else %otherwise
            adjacent=yI(2)-yI(1); %calculate lenght of the opposite side
            opposite=xI(2)-xI(1); %calculate length of the adjacent side
            angle=atand(opposite/adjacent); %calcualte the angle
            original=imrotate(original,-angle,'bilinear','crop'); %rotate the image
        end
        axes(handles.OriginalImageCroppedImage); %set to show on axes1
        imshow(original); %display image
    end
end



% --- Executes on button press in DetectViable.
function DetectViable_Callback(hObject, eventdata, handles)
% hObject    handle to DetectViable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global original DETECTViable_check DIL CROPIMAGE_check removebkgnd ...
    DETECTViableTHRESH_check FolderSelect_check CreateStruct FontSIZE ...
    image_crop

if FolderSelect_check==0 ||DETECTViableTHRESH_check==0
    if FolderSelect_check==0;
        msgbox(['\fontsize{' ,num2str(FontSIZE),'}Please select a folder.'], 'Warning','warn',CreateStruct);
    elseif FolderSelect_check==1;
        if DETECTViableTHRESH_check==0;
            msgbox(['\fontsize{' ,num2str(FontSIZE),'}Please detect the threshold for the border of the viable sample.'], 'Warning','warn',CreateStruct);
        end
    end
else
    DETECTViable_check=1;
    if CROPIMAGE_check==0
        redChannel=original(:, :, 1);
        greenChannel=original(:, :, 2);
        blueChannel=original(:, :, 3);
    elseif CROPIMAGE_check==1
        redChannel=image_crop(:, :, 1);
        greenChannel=image_crop(:, :, 2);
        blueChannel=image_crop(:, :, 3);
    end
    redChannel(~DIL) = 255;
    greenChannel(~DIL) = 255;
    blueChannel(~DIL) = 255;
    removebkgnd = cat(3, redChannel, greenChannel, blueChannel);
    axes(handles.OriginalImageCroppedImage)
    imshow(removebkgnd)
end




% --- Executes on button press in DetectNecrosis.
function DetectNecrosis_Callback(hObject, eventdata, handles)
% hObject    handle to DetectNecrosis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global DETECTViable_check DIL DILn necfrac DetectNecrosis_check necroticpixels ...
    viablepixels FolderSelect_check DetectNecrosisTHRESH_check DETECTNecrosis_check ...
    DETECTViableTHRESH_check black2 CreateStruct FontSIZE

if FolderSelect_check==0 || DETECTViable_check==0 || DetectNecrosisTHRESH_check==0 || DETECTViableTHRESH_check==0
    if FolderSelect_check==0;
        msgbox(['\fontsize{' ,num2str(FontSIZE),'}Please select a folder.', 'Warning','warn'],CreateStruct);
    elseif FolderSelect_check==1;
        if DETECTViableTHRESH_check==0;
            msgbox(['\fontsize{' ,num2str(FontSIZE),'}Please detect threshold for the viable cells.'], 'Warning','warn',CreateStruct);
        elseif DETECTViableTHRESH_check==1;
            if DETECTViable_check==0;
                msgbox(['\fontsize{' ,num2str(FontSIZE),'}Please detect the border of the viable sample.'], 'Warning','warn',CreateStruct);
            elseif DETECTViable_check==1;
                DetectNecrosis_check=1;
                if DetectNecrosisTHRESH_check==0
                    msgbox(['\fontsize{' ,num2str(FontSIZE),'}Please detect the threshold for necrosis.'], 'Warning','warn',CreateStruct);
                end
            end
        end
    end
else
    DETECTNecrosis_check=1;
    necroticpixels=sum(sum(DILn==1));
    viablepixels=sum(sum(DIL==1));
    necfrac=necroticpixels/viablepixels;
    axes(handles.NecroticFractionAxes); %set to show on axes1
    imshow(black2);
    name=regexprep('-%','-',num2str(necfrac*100)); %set filename to be name
    hold on
    text(5,3,name,'FontWeight','bold','HorizontalAlignment','left','FontSize',11,'Interpreter','none','Color','w'); %display the filename
    hold off
end



% --- Executes on slider movement.
function Viable_Callback(hObject, eventdata, handles)
% hObject    handle to Viable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global image_crop original CROPIMAGE_check BWopen image_blue_binary black ...
    DIL DETECTViableTHRESH_check FolderSelect_check CreateStruct FontSIZE
if FolderSelect_check==0;
    msgbox(['\fontsize{' ,num2str(FontSIZE),'}Please select a folder.'], 'Warning','warn',CreateStruct);
elseif FolderSelect_check==1;
    axes(handles.NADHaxes);
    imshow(black)
    axes(handles.ViableBorder);
    imshow(black)
    T=get(handles.Viable,'Value');
    if CROPIMAGE_check==0
        blue=original(:,:,1);
    elseif CROPIMAGE_check==1
        blue=image_crop(:,:,1);
    end
    E = entropyfilt(blue);
    Eim = rescale(E);
    image_blue_binary=imbinarize(Eim,T/255);
    BWopen = bwareaopen(image_blue_binary,5000);
    SE = strel(true(15));
    ER=imerode(BWopen,SE);
    f_ER=imfill(ER,'holes');
    DIL=imdilate(f_ER,SE);
    axes(handles.NADHaxes);
    imshow(DIL)
    
    [B_tissue,L_tissue] = bwboundaries(DIL,'noholes'); %detect binary boundaries
    stats = regionprops(L_tissue,'all'); %measure properties of image regions
    allareas=[stats.Area]; %pull out the areas as a matrix
    [M,I]=max(allareas); %find the maximum spot
    axes(handles.ViableBorder)
    if CROPIMAGE_check==1
        imshow(image_crop);
    elseif CROPIMAGE_check==0
        imshow(original);
    end
    hold on
    for k=1:length(B_tissue)
        boundary_tissue=B_tissue{k};
        plot(boundary_tissue(:,2),boundary_tissue(:,1),'g','linewidth',1);
    end
    drawnow expose
    DETECTViableTHRESH_check=1;
    hold off
end

% --- Executes during object creation, after setting all properties.
function Viable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Viable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function NecroticThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to NecroticThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global original DILn removebkgnd DetectNecrosisTHRESH_check DETECTViable_check ...
    DETECTViableTHRESH_check FolderSelect_check CreateStruct FontSIZE CROPIMAGE_check ...
    image_crop

if FolderSelect_check==0 || DETECTViable_check==0 ||DETECTViableTHRESH_check==0
    if FolderSelect_check==0;
        msgbox(['\fontsize{' ,num2str(FontSIZE),'}Please select a folder.'], 'Warning','warn',CreateStruct);
    elseif FolderSelect_check==1;
        if DETECTViableTHRESH_check==0;
            msgbox(['\fontsize{' ,num2str(FontSIZE),'}Please detect threshold for the viable cells.'], 'Warning','warn',CreateStruct);
        elseif DETECTViableTHRESH_check==1;
            if DETECTViable_check==0;
                msgbox(['\fontsize{' ,num2str(FontSIZE),'}Please detect the border of the viable sample.'], 'Warning','warn',CreateStruct);
            end
        end
    end
else
    DetectNecrosisTHRESH_check=1;
    T2=get(handles.NecroticThreshold,'Value')
    blueN=removebkgnd(:,:,1);
    axes(handles.OriginalImageCroppedImage);
    imshow(removebkgnd);
    EN = entropyfilt(blueN);
    EimN = rescale(EN);
    EimN(EimN==0)=255;
    image_blue_binaryN=imbinarize(EimN,T2/255);
    BWopenN = bwareaopen(~image_blue_binaryN,2000);
    SE = strel(true(15));
    ER=imerode(BWopenN,SE);
    f_ER_N=imfill(ER,'holes');
    DILn=imdilate(f_ER_N,SE);
    axes(handles.NecroticFraction)
    imshow(DILn);
    [B_necrotic,L_necrotic] = bwboundaries(DILn,'noholes'); %detect binary boundaries
    stats = regionprops(L_necrotic,'all'); %measure properties of image regions
    allareas=[stats.Area]; %pull out the areas as a matrix
    [M,I]=max(allareas); %find the maximum spot
    axes(handles.NecroticBorder)
    if CROPIMAGE_check==1
    imshow(image_crop);
    elseif CROPIMAGE_check==0
    imshow(original);
    end
    hold on
    for k=1:length(B_necrotic)
        boundary_necrotic=B_necrotic{k};
        plot(boundary_necrotic(:,2),boundary_necrotic(:,1),'r','linewidth',1);
    end
    drawnow expose
    hold off
end


% --- Executes during object creation, after setting all properties.
function NecroticThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NecroticThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in DisplayOriginal.
function DisplayOriginal_Callback(hObject, eventdata, handles)
% hObject    handle to DisplayOriginal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global original original_noedits CROPIMAGE_check black  black2 ...
    SAVE_check NextImage_Check FolderSelect_check CreateStruct FontSIZE
if FolderSelect_check==0;
    msgbox(['\fontsize{' ,num2str(FontSIZE),'}Please select a folder.'], 'Warning','warn',CreateStruct);
elseif FolderSelect_check==1;
    original=original_noedits;
    CROPIMAGE_check=0;
    SAVE_check=0;
    NextImage_Check=0;
    axes(handles.NADHaxes); %set to show on axes1
    imshow(black); %display image
    axes(handles.NecroticFraction); %set to show on axes1
    imshow(black); %display image
    axes(handles.ViableBorder); %set to show on axes1
    imshow(black); %display image
    axes(handles.OriginalImageCroppedImage); %set to show on axes1
    imshow(original); %display image
end


% --- Executes on button press in Save_Data.
function Save_Data_Callback(hObject, eventdata, handles)
% hObject    handle to Save_Data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global  SAVE_check  NextImage_Check filename  ...
    necfrac necroticpixels viablepixels DATA i FolderSelect_check ...
    DETECTViable_check DetectNecrosisTHRESH_check DETECTNecrosis_check ...
    DETECTViableTHRESH_check CreateStruct FontSIZE
if FolderSelect_check==0 || DETECTViable_check==0 || DetectNecrosisTHRESH_check==0 || DETECTNecrosis_check==0 
    if FolderSelect_check==0;
        msgbox(['\fontsize{' ,num2str(FontSIZE),'}Please select a folder.'], 'Warning','warn',CreateStruct);
    elseif FolderSelect_check==1;
        if DETECTViableTHRESH_check==0;
            msgbox(['\fontsize{' ,num2str(FontSIZE),'}Please detect threshold for the viable cells.'], 'Warning','warn',CreateStruct);
        elseif DETECTViableTHRESH_check==1;
            if DETECTViable_check==0;
                msgbox(['\fontsize{' ,num2str(FontSIZE),'}Please detect the border of the viable sample.'], 'Warning','warn',CreateStruct);
            elseif DETECTViable_check==1;
                if DetectNecrosisTHRESH_check==0
                    msgbox(['\fontsize{' ,num2str(FontSIZE),'}Please detect the threshold for necrosis.'], 'Warning','warn',CreateStruct);
                elseif DetectNecrosisTHRESH_check==1
                    if DETECTNecrosis_check==0
                        msgbox(['\fontsize{' ,num2str(FontSIZE),'}Please detect the necrotic fraction.'], 'Warning','warn',CreateStruct);
                    end
                end
            end
        end
    end
else
    SAVE_check=1;
    DATA{i,1}=filename;
    DATA{i,2}=necfrac;
    DATA{i,3}=necroticpixels;
    DATA{i,4}=viablepixels;
    NextImage_Check=0;
end


% --- Executes on button press in NextImage.
function NextImage_Callback(hObject, eventdata, handles)
% hObject    handle to NextImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global  pathname  original_noedits  filename ...
    CROPIMAGEcheck SAVE_check original NextImage_Check i FolderSelect_check ...
    DETECTViable_check DetectNecrosisTHRESH_check DETECTNecrosis_check filename ...
    files DETECTViableTHRESH_check CROPIMAGE_check CreateStruct FontSIZE black black2
if FolderSelect_check==0 || DETECTViable_check==0 || DetectNecrosisTHRESH_check==0 || DETECTNecrosis_check==0
        SAVE_check==0 
    if FolderSelect_check==0;
        msgbox(['\fontsize{' ,num2str(FontSIZE),'}Please select a folder.'], 'Warning','warn',CreateStruct);
    elseif FolderSelect_check==1;
        if DETECTViableTHRESH_check==0;
            msgbox(['\fontsize{' ,num2str(FontSIZE),'}Please detect threshold for the viable cells.'], 'Warning','warn',CreateStruct);
        elseif DETECTViableTHRESH_check==1;
            if DETECTViable_check==0;
                msgbox(['\fontsize{' ,num2str(FontSIZE),'}Please detect the border of the viable sample.'], 'Warning','warn',CreateStruct);
            elseif DETECTViable_check==1;
                if DetectNecrosisTHRESH_check==0
                    msgbox(['\fontsize{' ,num2str(FontSIZE),'}Please detect the threshold for necrosis.'], 'Warning','warn',CreateStruct);
                elseif DetectNecrosisTHRESH_check==1
                    if DETECTNecrosis_check==0
                        msgbox(['\fontsize{' ,num2str(FontSIZE),'}Please detect the necrotic fraction.'], 'Warning','warn',CreateStruct);
                    elseif DETECTNecrosis_check==1
                        if  SAVE_check==0
                            msgbox(['\fontsize{' ,num2str(FontSIZE),'}Please save the results.'], 'Warning','warn',CreateStruct);
                        end
                    end
                end
            end
        end
    end
else
    path ='Y:\data\Erika\NADH Save Space\';
    F_NADH=getframe(handles.ViableBorder);
    F_Necrotic=getframe(handles.NecroticBorder);
    Tissue=frame2im(F_NADH);
    Necrosis=frame2im(F_Necrotic);
    filename_notiff=filename(1:end-5);
    imwrite(Tissue,[path filename_notiff '_Sample.tiff' ]);
    imwrite(Necrosis,[path filename_notiff '_Necrosis.tiff' ]);
    if i>=length(files)
        msgbox(['\fontsize{' ,num2str(FontSIZE),'}End of Folder. Export results.'], 'Warning','warn',CreateStruct);                   
    else
        
    CROPIMAGEcheck=0;
    NextImage_Check=1;
    i=i+1;
    filename=files(i).name;
    original=imread([pathname '\' filename]);
    axes(handles.OriginalImageCroppedImage); %set to show on axes1
    
    imshow(original); %display image
    axes(handles.NADHaxes); %set to show on axes1
    
    imshow(black); %display image
    axes(handles.NecroticFraction); %set to show on axes1
    
    imshow(black); %display image
    axes(handles.ViableBorder); %set to show on axes1
    
    imshow(black); %display image
    axes(handles.NecroticBorder);
    
    imshow(black)
    axes(handles.NecroticFractionAxes);
    imshow(black2);
    original_noedits=original;
    SAVE_check=0;
    CROPIMAGE_check=0;
    DETECTViable_check=0;
    DETECTNecrosis_check=0;
    DETECTViableTHRESH_check=0;
    DetectNecrosisTHRESH_check=0;
    end
end




% --- Executes on button press in ExportData.
function ExportData_Callback(hObject, eventdata, handles)
% hObject    handle to ExportData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global DATA FolderSelect_check DETECTViable_check DetectNecrosisTHRESH_check ...
    DETECTNecrosis_check files DETECTViableTHRESH_check SAVE_check i ...
    CROPIMAGE_check CreateStruct FontSIZE
if FolderSelect_check==0 || DETECTViable_check==0 || DetectNecrosisTHRESH_check==0 || DETECTNecrosis_check==0 || ...
        SAVE_check==0 || i<length(files)
    if FolderSelect_check==0
        msgbox(['\fontsize{' ,num2str(FontSIZE),'}Please select a folder.'], 'Warning','warn',CreateStruct);
    elseif FolderSelect_check==1;
        if DETECTViable_check==0;
            msgbox(['\fontsize{' ,num2str(FontSIZE),'}Please detect threshold for the viable cells.'], 'Warning','warn',CreateStruct);
        elseif DETECTViableTHRESH_check==1;
            if DETECTViable_check==0;
                msgbox(['\fontsize{' ,num2str(FontSIZE),'}Please detect the border of the viable sample.'], 'Warning','warn',CreateStruct);
            elseif DETECTViable_check==1;
                if DetectNecrosisTHRESH_check==0
                    msgbox(['\fontsize{' ,num2str(FontSIZE),'}Please detect the threshold for necrosis.'], 'Warning','warn',CreateStruct);
                elseif DetectNecrosisTHRESH_check==1
                    if DETECTNecrosis_check==0
                        msgbox(['\fontsize{' ,num2str(FontSIZE),'}Please detect the necrotic fraction.'], 'Warning','warn',CreateStruct);
                    elseif DETECTNecrosis_check==1
                        if  SAVE_check==0
                            msgbox(['\fontsize{' ,num2str(FontSIZE),'}Please save the results.'], 'Warning','warn',CreateStruct);
                        else
                            if i<length(files)
                                msgbox(['\fontsize{' ,num2str(FontSIZE),'}Please finish all images in the folder.'], 'Warning','warn',CreateStruct);
                            end
                        end
                    end
                end
            end
        end
    end
else
    prev_DATA=readcell('NecroticFractions.xls')
    DATA
    save_DATA=[prev_DATA; DATA];
    save_DATA
    writecell(save_DATA,'NecroticFractions.xls')
    SAVE_check=0;
    CROPIMAGE_check=0;
    DETECTViable_check=0;
    DETECTNecrosis_check=0;
    SAVE_check=0;
    FolderSelect_check=0;
    DETECTViableTHRESH_check=0;
    DetectNecrosisTHRESH_check=0;
end



% --- Executes on button press in SelectFolder.
function SelectFolder_Callback(hObject, eventdata, handles)
% hObject    handle to SelectFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ID pathname  stackL stackFL stackBoth  original black original_noedits  ...
    CROPIMAGE_check SAVE_check i files filename FolderSelect_check CreateStruct FontSIZE ...
    black black2 DATA
DATA={};
CROPIMAGE_check=0;
SAVE_check=0;
FolderSelect_check=1;
i=1;
axes(handles.OriginalImageCroppedImage); %set to show on axes1
imshow(black); %display image
axes(handles.NADHaxes); %set to show on axes1
imshow(black); %display image
axes(handles.NecroticFraction); %set to show on axes1
imshow(black); %display image
axes(handles.ViableBorder); %set to show on axes1
imshow(black); %display image
axes(handles.NecroticBorder); %set to show on axes1
imshow(black); %display image
axes(handles.FolderName); %set to show on axes1
imshow(black2);
pathname = uigetdir('Y:\data\Erika\Pathology Images Tiff\','Select Folder');
files=dir(fullfile(pathname, '*.tiff'));
slashes=find(pathname=='\');
ID=pathname(slashes(end)+1:end);
axes(handles.FolderName); %set to show on axes1
imshow(black2); %display image
name=regexprep('Folder: -','-',ID); %set filename to be name
hold on
text(5,3,name,'FontWeight','bold','HorizontalAlignment','left','FontSize',11,'Interpreter','none','Color','w'); %display the filename
hold off
filename=files(i).name;
original=imread([pathname '\' filename]);
axes(handles.OriginalImageCroppedImage); %set to show on axes1
imshow(original); %display image
original_noedits=original;
stackL=[];
stackFL=[];
stackBoth=[];
