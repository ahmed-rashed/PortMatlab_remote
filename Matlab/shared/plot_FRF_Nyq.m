function curve_handle=plot_FRF_Nyq(H_vec, ...
                         H_label,DispMagLines,varargin)  %Optional arguments

H_real_multiplier='';
if nargin<2
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
H_Latex_subtitle=[H_Latex_subtitle,'\left(f \right)'];

curve_handle=plot(real(H_vec),imag(H_vec),varargin{:});
if nargin>2
    if DispMagLines
        n_f=length(H_vec);
        n_MagLines=DispMagLines;
        delta_temp=floor(n_f/n_MagLines);
        iidx=1:delta_temp:n_f;

        set(gcf,'NextPlot','add');
        set(gca,'NextPlot','add');
        setappdata(gca,'PlotHoldStyle',false);

        if n_f==size(H_vec,1)
            H_temp_vec=H_vec(iidx,:);
        else     %Assuming that n_f==size(H_vec,2)
            H_temp_vec=H_vec(:,iidx);
        end

        handle=quiver(0*H_temp_vec,0*H_temp_vec,real(H_temp_vec),imag(H_temp_vec),0,'ShowArrowHead','off','LineWidth',get(curve_handle,'LineWidth')/3);
        set(get(get(handle,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
        
        handle2=plot(real(H_temp_vec),imag(H_temp_vec),'.','MarkerSize',15,'MarkerEdgeColor',get(handle,'Color'));
        set(get(get(handle2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
    end
end

vax=axis;
DataAspectRatio=get(gca,'DataAspectRatio');
set(gca,'DataAspectRatio',[DataAspectRatio(1),DataAspectRatio(1),DataAspectRatio(3)]);
axis(vax);

xlabel(['$\Re\left(',H_Latex_subtitle,'\right)',H_real_multiplier,'$'], 'interpreter', 'latex')
ylabel(['$\Im\left(',H_Latex_subtitle,'\right)',H_real_multiplier,'$'], 'interpreter', 'latex')
grid on