function [] = sort_licks_gui_run(varargin)
if (nargin == 0)
    flag = 'start';
else
    flag = varargin{1};
end

switch flag
    case 'start'
        local_start;
    case 'load'
        local_load;
end

end

function local_start()
% open gui
sort_licks_gui;

% intialize strings
set(findobj('tag','keep_button'),'string','KEEP')
set(findobj('tag','reject_button'),'string','REJECT')
set(findobj('tag','load_data_button'),'string','Load Data')

% set callbacks
set(findobj('tag','load_data_button'),'callback',['sort_licks_gui_run(''load'')'])
end

function local_load()
global licks curLick curChan
[fname,~] = uigetfile('*.mat');
licks = load(fname);
licks = licks.licks;
axes_name = findobj('tag','data_axes');
fig = findobj('tag','sort_licks_gui');
set(fig,'CurrentAxes',axes_name)
curChan = 1;
curLick = 1;
plot(licks{curChan}(curLick).data)
end