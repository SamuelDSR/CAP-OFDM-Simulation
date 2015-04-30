function [xdata,ydata] = get_xy_from_figures()
% Get a handle to the current figure:
h = gcf; %current figure handle   
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes   

%objTypes = get(dataObjs, 'Type');  %type of low-level graphics object   
%NOTE : Different objects like 'Line' and 'Surface' will store data differently. Based on the 'Type', you can search the documentation for how each type stores its data.

xdata = get(dataObjs, 'XData');  %data from low-level grahics objects

ydata = get(dataObjs, 'YData');

%zdata = get(dataObjs, 'ZData');
