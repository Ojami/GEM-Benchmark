function methodChoice = methodSelection(handles)

methodChoice = []; cnt = 1;
checkVal = get(handles.pFBAcheck,'Value');
if checkVal
    methodChoice{cnt} = get(handles.pFBAcheck,'String');
    cnt = cnt + 1;
end
clear checkVal
checkVal = get(handles.TRFBAcheck,'Value');
if checkVal
    methodChoice{cnt} = get(handles.TRFBAcheck,'String');
    cnt = cnt + 1;
end
clear checkVal
checkVal = get(handles.PRIMEcheck,'Value');
if checkVal
    methodChoice{cnt} = get(handles.PRIMEcheck,'String');
    cnt = cnt + 1;
end
clear checkVal
checkVal = get(handles.GIMMEcheck,'Value');
if checkVal
    methodChoice{cnt} = get(handles.GIMMEcheck,'String');
    cnt = cnt + 1;
end
clear checkVal
checkVal = get(handles.iMATcheck,'Value');
if checkVal
    methodChoice{cnt} = get(handles.iMATcheck,'String');
    cnt = cnt + 1;
end
clear checkVal
checkVal = get(handles.INITcheck,'Value');
if checkVal
    methodChoice{cnt} = get(handles.INITcheck,'String');
    cnt = cnt + 1;
end
clear checkVal
checkVal = get(handles.mCADREcheck,'Value');
if checkVal
    methodChoice{cnt} = get(handles.mCADREcheck,'String');
    cnt = cnt + 1;
end
clear checkVal
checkVal = get(handles.FASTCOREcheck,'Value');
if checkVal
    methodChoice{cnt} = get(handles.FASTCOREcheck,'String');
    cnt = cnt + 1;
end
clear checkVal
checkVal = get(handles.FASTCORMICScheck,'Value');
if checkVal
    methodChoice{cnt} = get(handles.FASTCORMICScheck,'String');
    cnt = cnt + 1;
end
clear checkVal
checkVal = get(handles.CORDAcheck,'Value');
if checkVal
    methodChoice{cnt} = get(handles.CORDAcheck,'String');
    cnt = cnt + 1;
end
clear checkVal
if isempty(methodChoice)
    uiwait(errordlg('You should select at least one algorithm!','Select Algorithm(s)'));
    return
end