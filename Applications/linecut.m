function [cuthandles] = linecut(plothandle)
%LINECUT View linecuts along the X and Y dimension of a 3D plot
%   LINECUT(H) places a crosshair on the 3D plot with handle H. H can be a
%   handle to a pcolor, surf, etc. plot or a non-RGB image. The crosshair
%   can be moved with either the mouse or the arrow keys. Pressing enter
%   removes the crosshair and stops linecut.
%   
%   CH = LINECUT(H) returns a two element vector of handles to the X and Y
%   linecut objects respectively
%
%   Note that surface plots must have X and Y data that are either vectors
%   or meshgrid formatted matrices.

%   Patrick Maher
%   v. 1.0 9/7/2014

% https://fr.mathworks.com/matlabcentral/fileexchange/47766-linecut

 % test if this is an object with XData
   has_xdata = false
   while ~has_xdata
    try
      get(plothandle,'Xdata');
      has_xdata = true;
      plothandle = plothandle(1);
    catch
      plothandle = get(plothandle, 'Children');
    end
  end
  hand.mainplot=plothandle;
  hand.mainax=get(hand.mainplot,'Parent');
  hand.mainax_title = get(hand.mainax,'Title');
  hand.mainax_title = get(hand.mainax_title, 'String');
  hand.mainax_xlabel = get(hand.mainax,'XLabel');
  hand.mainax_xlabel = get(hand.mainax_xlabel, 'String');
  hand.mainax_ylabel = get(hand.mainax,'YLabel');
  hand.mainax_ylabel = get(hand.mainax_ylabel, 'String');
  hand.mainfig=get(hand.mainax,'Parent');
  hand.Xdata=get(hand.mainplot,'Xdata');
  hand.Ydata=get(hand.mainplot,'Ydata');
  %make sure X and Y data either vectors or product of meshgrid
  if ~isvector(hand.Xdata) || ~isvector(hand.Ydata)
      if isequal(hand.Xdata,repmat(hand.Xdata(1,:),size(hand.Xdata,1),1)) && isequal(hand.Ydata,repmat(hand.Ydata(:,1),1,size(hand.Ydata,2)))
          hand.Xdata=hand.Xdata(1,:);
          hand.Ydata=hand.Ydata(:,1);
      else
          error('X and Y data must square');
      end
  end

  if strcmp(get(hand.mainplot,'Type'),'image') %if object is an image
      hand.Zdata=get(hand.mainplot,'Cdata'); 
      if length(size(hand.Zdata))>2
          error('RGB data not supported');
      end
      hand.plottype='image';
      %make sure X,Y vectors are same size as image data
      if length(hand.Xdata)==1
          hand.Xdata=hand.Xdata:(hand.Xdata+size(hand.Zdata,2)-1);
      else
          hand.Xdata=linspace(hand.Xdata(1),hand.Xdata(end),size(hand.Zdata,2));
      end
      if length(hand.Ydata)==1
          hand.Ydata=hand.Ydata:(hand.Ydata+size(hand.Zdata,1)-1);
      else
          hand.Ydata=linspace(hand.Ydata(1),hand.Ydata(end),size(hand.Zdata,1));
      end
  elseif all(all(get(hand.mainplot,'Zdata')==0)) %if object is pcolor
      hand.Zdata=get(hand.mainplot,'Cdata');
      hand.plottype='pcolor';
  else
      hand.Zdata=get(hand.mainplot,'Zdata');
      hand.plottype='surf';
  end

  hand.xindex=1;hand.yindex=1; %initial indices
  [hand.xcutplot, hand.ycutplot]=deal([]); %declare variables
  %draw crosshairs
  if any(strcmp(hand.plottype,{'image','pcolor'}))
      hand.lx=line(xlim(hand.mainax),[hand.Ydata(1) hand.Ydata(1)]);
      hand.ly=line([hand.Xdata(1) hand.Xdata(1)],ylim(hand.mainax));
      set([hand.lx hand.ly],'LineWidth',2);
  else
      hand.lx=patch([hand.Xdata(1) hand.Xdata(1) hand.Xdata(end) hand.Xdata(end)],[hand.Ydata(1) hand.Ydata(1) hand.Ydata(1) hand.Ydata(1)],[zlim(hand.mainax) fliplr(zlim(hand.mainax))],'b');
      hand.ly=patch([hand.Xdata(1) hand.Xdata(1) hand.Xdata(1) hand.Xdata(1)],[hand.Ydata(1) hand.Ydata(1) hand.Ydata(end) hand.Ydata(end)],[zlim(hand.mainax) fliplr(zlim(hand.mainax))],'b');
      hand.lintersect=line([hand.Xdata(1) hand.Xdata(1)], [hand.Ydata(1) hand.Ydata(1)], zlim(hand.mainax),'Color','k');
      set(hand.lintersect,'Parent',hand.mainax,'HitTest','off');
      set([hand.lx hand.ly],'FaceAlpha',0.2);
  end
  %make sure clicking passes through crosshairs
  set(hand.lx,'Parent',hand.mainax,'HitTest','off');
  set(hand.ly,'Parent',hand.mainax,'HitTest','off');

  set(hand.mainfig,'WindowButtonDownFcn',@MouseClick, ...
    'KeyPressFcn',@KeyPress, 'CloseRequestFcn', @CloseRequestFcn);
  guidata(hand.mainplot,hand);
  cuthandles=updatecuts(hand.mainfig);
end % main linecut



function MouseClick(src,EventData)
  hand=guidata(src);
  point=get(hand.mainax,'CurrentPoint');
  %find indices closest to click
  [~,hand.xindex]=min(abs(hand.Xdata-point(1,1)));
  [~,hand.yindex]=min(abs(hand.Ydata-point(1,2)));
  guidata(hand.mainfig, hand);
  updatecuts(hand.mainfig);
end % MouseClick

function KeyPress(src,EventData)
  hand=guidata(src);
  hand=guidata(hand.mainfig);
  key=EventData.Key;
  xindinc=1;yindinc=1;
  %make sure arrow key direction corresponds to crosshair motion
  if all(diff(hand.Xdata)<0), xindinc=-xindinc; end
  if all(diff(hand.Ydata)<0), yindinc=-yindinc; end
  if strcmp(get(hand.mainax,'Xdir'),'reverse'), xindinc=-xindinc; end
  if strcmp(get(hand.mainax,'Ydir'),'reverse'), yindinc=-yindinc; end

  switch key
      case 'rightarrow'
          hand.xindex = hand.xindex+xindinc;
      case 'leftarrow'
          hand.xindex = hand.xindex-xindinc;
      case 'uparrow'
          hand.yindex = hand.yindex+yindinc;
      case 'downarrow'
          hand.yindex = hand.yindex-yindinc;
      case 'return'
          delete(hand.lx);
          delete(hand.ly);
          if isfield(hand,'lintersect'), delete(hand.lintersect); end
          set(hand.mainfig,'WindowButtonDownFcn',{},'KeyPressFcn',{});       
          set(hand.xcutfig,'KeyPressFcn',{});
          set(hand.ycutfig,'KeyPressFcn',{});
          return;        
  end
  %make sure indices aren't out of bounds
  hand.xindex=max(hand.xindex,1); hand.xindex=min(hand.xindex,length(hand.Xdata));
  hand.yindex=max(hand.yindex,1); hand.yindex=min(hand.yindex,length(hand.Ydata));
  guidata(hand.mainfig,hand);
  updatecuts(hand.mainfig);
end % KeyPress


function cuthandles=updatecuts(mainfig)
  hand=guidata(mainfig);
  %update linecut plots
  if ishghandle(hand.xcutplot)
      set(hand.xcutplot,'YData', hand.Zdata(:,hand.xindex));
  else
      hand.xcutfig=figure;
      hand.xcutplot=plot(hand.Ydata, hand.Zdata(:,hand.xindex));
      xlabel(hand.mainax_ylabel);
      hand.xcutax=gca;
      set(hand.xcutfig,'KeyPressFcn',@KeyPress)
      guidata(hand.xcutfig,hand);
  end
  if ishghandle(hand.ycutplot)
      set(hand.ycutplot,'YData', hand.Zdata(hand.yindex,:));
  else
      hand.ycutfig=figure;
      hand.ycutplot=plot(hand.Xdata, hand.Zdata(hand.yindex,:));
      xlabel(hand.mainax_xlabel);
      hand.ycutax=gca;
      set(hand.ycutfig,'KeyPressFcn',@KeyPress)
      guidata(hand.ycutfig,hand);
  end
  %move crosshairs
  if any(strcmp(hand.plottype,{'image','pcolor'}))
      set(hand.lx,'YData',[hand.Ydata(hand.yindex) hand.Ydata(hand.yindex)]);
      set(hand.ly,'XData',[hand.Xdata(hand.xindex) hand.Xdata(hand.xindex)]);
  else
      set(hand.lx,'YData',hand.Ydata(hand.yindex)*ones(1,4));
      set(hand.ly,'XData',hand.Xdata(hand.xindex)*ones(1,4));
      set(hand.lintersect,'XData',[hand.Xdata(hand.xindex) hand.Xdata(hand.xindex)],'YData',[hand.Ydata(hand.yindex) hand.Ydata(hand.yindex)]);
  end

  set(get(hand.ycutax,'Title'),'String',[hand.mainax_title ' y=' num2str(hand.Ydata(hand.yindex)) ' (index #' num2str(hand.yindex) ')']);
  set(get(hand.xcutax,'Title'),'String',[hand.mainax_title ' x=' num2str(hand.Xdata(hand.xindex)) ' (index #' num2str(hand.xindex) ')']);
  guidata(hand.mainplot, hand);
  %return cut plot handles
  cuthandles=[hand.xcutplot hand.ycutplot];
end % updatecuts

function CloseRequestFcn(src,EventData)
  hand=guidata(src);
  hand=guidata(hand.mainfig);
  % close all
  delete(hand.lx);
  delete(hand.ly);
  delete(hand.xcutfig);
  delete(hand.ycutfig);
  delete(hand.mainfig);
end % CloseRequestFcn
