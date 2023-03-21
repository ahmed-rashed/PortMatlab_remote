function [ax_mag_h,ax_phase_h,curve_mag_h,curve_phase_h]= ...
         plot_FRF_mag_phase(f_vec,H_vec, ...
                            islin,ax_mag_h,ax_phase_h,f_label,H_label,DispMagLines,maxPhaseLag, ... %Optional arguments
                            varargin)
%To create a nonvalid graphics object (axis handle), use eigther:
%1) h=gobjects
%2) h=matlab.graphics.GraphicsPlaceholder

if nargin<3
    islin=true;
else
    if isempty(islin)
        islin=true;
    end
end

if nargin<4
    ax_mag_h=subplot(4,1,[1:3]);
else
    if isempty(ax_mag_h)
        ax_mag_h=subplot(4,1,[1:3]);
    end
end
    
if nargin<5
    ax_phase_h=subplot(4,1,4);
else
    if isempty(ax_phase_h)
        ax_phase_h=subplot(4,1,4);
    end
end

if nargin<6
    f_label='$f$ (Hz)';
else
    if isempty(f_label)
        f_label='$f$ (Hz)';
    end
end
indices=strfind(f_label,'$');
if length(indices)<2,error('f_label does not include LaTeX inline equation !!'),end

H_real_multiplier='';
if nargin<7
    H_Latex_subtitle='H';
else
    if isempty(H_label)
        H_Latex_subtitle='H';
    elseif iscellstr(H_label)
        if length(H_label)~=2
            error('If H_label is cell string, it must have two elements; one for H_Latex_subtitle and the other for H_real_multiplier')
        end
        H_Latex_subtitle=H_label{1};
        H_real_multiplier=H_label{2};
    else
        H_Latex_subtitle=H_label;
    end
end
H_Latex_subtitle=[H_Latex_subtitle,'\left(',f_label(indices(1)+1:indices(2)-1),'\right)'];

curve_mag_h=plot(ax_mag_h,f_vec,abs(H_vec),varargin{:});
n_f=length(f_vec);
        
if nargin>7
    if DispMagLines
        n_MagLines=DispMagLines;
        delta_temp=floor(n_f/n_MagLines);
        iidx=1:delta_temp:n_f;

        set(gcf,'NextPlot','add');
        set(gca,'NextPlot','add');
        setappdata(gca,'PlotHoldStyle',false);

        if n_f==size(H_vec,1)
            H_temp=H_vec(iidx,:);
        else     %Assuming that n_f==size(H_vec,2)
            H_temp=H_vec(:,iidx);
        end

        handle=stem(ax_mag_h,f_vec(iidx),abs(H_temp),'LineWidth',get(curve_mag_h,'LineWidth')/3,'Marker','none');
        set(get(get(handle,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
        
        handle2=plot(ax_mag_h,f_vec(iidx),abs(H_temp),'.','MarkerSize',15,'MarkerEdgeColor',get(handle,'Color'));
        set(get(get(handle2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
    end
end

if nargin<9
	maxPhaseLag=[];
end

if islin(1)
    set(ax_mag_h,'YScale','lin')
else
    set(ax_mag_h,'YScale','log')
end
if length(islin)>1
    if islin(2)
        set(ax_mag_h,'XScale','lin')
        if isgraphics(ax_phase_h),set(ax_phase_h,'XScale','lin'),end
    else
        set(ax_mag_h,'XScale','log')
        if isgraphics(ax_phase_h),set(ax_phase_h,'XScale','log'),end
    end
end

if isgraphics(ax_phase_h)
	set(ax_mag_h,'XTickLabel',[]);
else
	xlabel(ax_mag_h,f_label, 'interpreter', 'latex')
end

% if iscell(H_Latex_subtitle)
%     legStr=cell(size(H_Latex_subtitle));
%     for ii=1:length(H_Latex_subtitle)
%         legStr{ii}=['$\left|',H_Latex_subtitle{ii},'\right|$'];
%     end
%     legend(ax_mag_h,legStr, 'interpreter', 'latex')
% else
    ylabel(ax_mag_h,['$\left|',H_Latex_subtitle,'\right|',H_real_multiplier,'$'], 'interpreter', 'latex')
% end
set(ax_mag_h,'xGrid','on','YGrid','on');

%Plot the phase
if isgraphics(ax_phase_h)
    H_angle_vec=correctedPhase(H_vec,maxPhaseLag);
    curve_phase_h=plot(ax_phase_h,f_vec,H_angle_vec,varargin{:});
    if nargin>7
        if DispMagLines
            set(gcf,'NextPlot','add');
            set(gca,'NextPlot','add');
            setappdata(gca,'PlotHoldStyle',false);

            if n_f==size(H_vec,1)
                H_angle_vec_temp=H_angle_vec(iidx,:);
            else     %Assuming that n_f==size(H_vec,2)
                H_angle_vec_temp=H_angle_vec(:,iidx);
            end

            handle2=plot(ax_phase_h,f_vec(iidx),H_angle_vec_temp,'.','MarkerSize',15,'MarkerEdgeColor',get(handle,'Color'));
            set(get(get(handle2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
        end
    end

    xlabel(ax_phase_h,f_label, 'interpreter', 'latex')
%     if iscell(H_Latex_subtitle)
%         for ii=1:length(H_Latex_subtitle)
%             legStr{ii}=['$\angle ',H_Latex_subtitle{ii},'$'];
%         end
%         legend(ax_phase_h,legStr, 'interpreter', 'latex')
%     else
        ylabel(ax_phase_h,['$\angle \left(',H_Latex_subtitle,'\right)$'], 'interpreter', 'latex')
%     end

    yStep=pi;
    %yTicks=get(ax_phase_h,'YTick');
    %set(ax_phase_h,'YTick',yStep*[floor(yTicks(1)/yStep):ceil(yTicks(end)/yStep)]);
    if strcmp(get(ax_phase_h,'YLimMode'),'manual')
        ylimits=get(ax_phase_h,'ylim');
        yLim_min=min(yStep*floor(min(min(H_angle_vec))/yStep),ylimits(1));
        yLim_max=max(yStep*ceil(max(max(H_angle_vec))/yStep),ylimits(2));
    else
        yLim_min=yStep*floor(min(min(H_angle_vec))/yStep);
        yLim_max=yStep*ceil(max(max(H_angle_vec))/yStep);
    end
    N_ticks=round((yLim_max-yLim_min)/yStep)+1;
    if N_ticks>6
        yStep=yStep*round(N_ticks/6);
    end
    yTicks=[yLim_min:yStep:yLim_max];
    N_ticks=length(yTicks);
    set(ax_phase_h,'YTick',yTicks);
    yTickLabels=strings(1,N_ticks);
    for ii=1:N_ticks
        if yTicks(ii)==0
            yTickLabels(ii)='$0$';
        elseif yTicks(ii)==pi
            yTickLabels(ii)='$\pi$';
        elseif yTicks(ii)==-pi
            yTickLabels(ii)='$-\pi$';
        else
            yTickLabels(ii)="$"+(yTicks(ii)/pi)+'\pi$';
        end
    end
    ax_phase_h.YAxis.TickLabelInterpreter='latex';
    set(ax_phase_h,'YTickLabel',yTickLabels);
    if yLim_min~=yLim_max
        ylim(ax_phase_h,[yLim_min,yLim_max])
    end
    set(ax_phase_h,'xGrid','on','YGrid','on');
    set(ax_phase_h,'YMinorTick','on');
end
