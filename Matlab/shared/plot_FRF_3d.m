function plot_FRF_3d(f_vec,H_vec, ...
                     f_label,H_label,DispRealImag,DispMagLines) %Optional arguments

curve_h=plot3(f_vec,real(H_vec),imag(H_vec));
view(30,20);
grid on

DataAspectRatio=get(gca,'DataAspectRatio');
set(gca,'DataAspectRatio',[DataAspectRatio(1),DataAspectRatio(2),DataAspectRatio(2)]);

if nargin<3
    f_label='$f$ (Hz)';
else
    if isempty(f_label)
        f_label='$f$ (Hz)';
    end
end
indices=strfind(f_label,'$');
if length(indices)<2,error('f_label does not include LaTeX inline equation !!'),end

H_real_multiplier='';
if nargin<4
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

xlabel(f_label,'interpreter', 'latex');
ylabel(['$\Re\left(',H_Latex_subtitle,'\right)',H_real_multiplier,'$'], 'interpreter', 'latex')
zlabel(['$\Im\left(',H_Latex_subtitle,'\right)',H_real_multiplier,'$'], 'interpreter', 'latex')

if nargin>4
    nexta = get(gca,'NextPlot');
    nextf = get(gcf,'NextPlot');

    if DispRealImag
        set(gca,'YLimMode','manual');
        set(gca,'ZLimMode','manual');
        ylims=ylim;
        zlims=zlim;

        set(gcf,'NextPlot','add');
        set(gca,'NextPlot','add');

        %Nyquist
        hLine=line(f_vec*0,real(H_vec),imag(H_vec),'Color',get(curve_h,'Color'),'LineStyle',get(curve_h,'LineStyle'),'LineWidth',get(curve_h,'LineWidth')/3);
        set(get(get(hLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend

        %Imag
        hLine=line(f_vec,max(ylims)*ones(size(f_vec)),imag(H_vec),'Color',get(curve_h,'Color'),'LineStyle',get(curve_h,'LineStyle'),'LineWidth',get(curve_h,'LineWidth')/3);
        set(get(get(hLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend

        %Real
        hLine=line(f_vec,real(H_vec),min(zlims)*ones(size(f_vec)),'Color',get(curve_h,'Color'),'LineStyle',get(curve_h,'LineStyle'),'LineWidth',get(curve_h,'LineWidth')/3);
        set(get(get(hLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
    end
    
    if nargin>5
        if DispMagLines
            n_f=length(f_vec);
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
            
            handle=quiver3(f_vec(iidx),0*f_vec(iidx),0*f_vec(iidx),0*f_vec(iidx),real(H_temp),imag(H_temp),0,'ShowArrowHead','off','LineWidth',get(curve_h,'LineWidth')/3);
            set(get(get(handle,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
            
            handle2=plot3(f_vec(iidx),real(H_temp),imag(H_temp),'.','MarkerSize',15,'MarkerEdgeColor',get(handle,'Color'));
            set(get(get(handle2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend

        end
    end

    set(gcf,'NextPlot',nexta);
    set(gca,'NextPlot',nextf);
end