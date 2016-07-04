function datetickzoom(varargin)
%DATETICKZOOM Date formatted tick labels, automatically updated when zoomed or panned.
%   Arguments are completely identical to does of DATETICK. 
%
%   See also datetick, datestr, datenum
%
% Works in matlab R2014b onwards. 
%
% EXAMPLE:
%   plot(now+(1:400),cos((1:400)*.2));
%   datetickzoom keeplimits
%
%
% Aslak Grinsted 2015. 

callbackdata=varargin;
%we need to ensure that a handle is passed to callbacks, so add gca if it
%is not there
if isempty(callbackdata) || ~all(isgraphics(callbackdata{1}))
    callbackdata={gca callbackdata{:}}; 
end
datetick(callbackdata{:})

%For callbacks add parameter keeplimits if it is not already there
ix=cellfun(@(x)strcmpi(x,'keeplimits'),callbackdata);
callbackdata(ix)=[];
callbackdata{end+1}='keeplimits'; %for callbacks use this always

addlistener(gcf,'SizeChanged',@(h,event,varargin)callback(h,event,callbackdata{:}));
addlistener(gca,'XLim','PostSet',@(h,event,varargin)callback(h,event,callbackdata{:}));
addlistener(gca,'YLim','PostSet',@(h,event,varargin)callback(h,event,callbackdata{:}));
addlistener(gca,'Position','PostSet',@(h,event,varargin)callback(h,event,callbackdata{:}));

function callback(h,event,varargin)
datetick(varargin{:})

