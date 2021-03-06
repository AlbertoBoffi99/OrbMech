% PORK-CHOP PLOTS FILE
% In this file pork-chop plots are plotted. If this is the first iteration
% graphs are created, otherwise they are just imported in Matlab figure
% format

%-------------------------------------------------------------------------%

if options.pcpfirsttime

%% PRELIMINARY ANALYSIS OF FIRST LAMBERT ARC

temp.n = 1000;
temp.nplanes = 5;

dep_window = linspace(departure.date_min, arrival.date_max, temp.n);
fb_window = linspace(departure.date_min, arrival.date_max, temp.n);

temp.Dvmin = astroConstants(5); 
r = NaN;
c = NaN;                                                  
% select an interval to cut the contour plot
planes = 1 : 1 : temp.nplanes;                                    

for i = 1 : temp.n
    % for every Jupiter's position at departure
    kep1 = uplanet(dep_window(i), 5);
    [r1(:,i), v1(:,i)] = kep2car(kep1(1),kep1(2),kep1(3),kep1(4),kep1(5),kep1(6),astro.muS);

    for j = 1 : temp.n
        %for every Earth's position at flyby
        kep2 = uplanet(fb_window(j), 3);
        [r2(:,j), v2(:,j)] = kep2car(kep2(1),kep2(2),kep2(3),kep2(4),kep2(5),kep2(6),astro.muS);
        
        T11(i,j) = (fb_window(j)-dep_window(i))*24*60*60;   
        %flyby corresponding to the i-departure date and j-flyby date
        [aT1(i,j),pT,E_T,ERROR,Vi1(i,j,:),Vf1(i,j,:),TPAR,THETA] = ...
            lambertMR(r1(:,i),r2(:,j),T11(i,j),astro.muS,options.orbitType,options.Nrev,options.Ncase,options.LambOptions);
        Dv1x = Vi1(i,j,1)-v1(1,i);
        Dv1y = Vi1(i,j,2)-v1(2,i);
        Dv1z = Vi1(i,j,3)-v1(3,i);
        Dv1 = norm([Dv1x Dv1y Dv1z]);  
        Dvtot(i,j) = Dv1;
        
        if Dvtot(i,j) < temp.Dvmin 
            temp.Dvmin = Dvtot(i,j);
            r = i;
            c = j;
        end
    end
    x(i) = datenum(mjd20002date(dep_window(i)));  
    y(i) = datenum(mjd20002date(fb_window(i))); 
end 

% dates obtained from overall problem
l_day = mjd20002date(results.dates(1));
fb_day = mjd20002date(results.dates(2));

pcpL1 = figure();
pcpL1.WindowState = 'maximized';
hold on
[C, h] = contour(x, y, Dvtot', floor(temp.Dvmin) + (planes));
plot(x(r), y(c), '*', 'LineWidth',12);
caxis([temp.Dvmin temp.Dvmin + planes(length(planes))]);                                     
Color = colorbar;  
colormap('Turbo');
clabel(C,h, floor(temp.Dvmin) + (planes));                         
datetick( 'x', 'yyyy mmm dd', 'keeplimits', 'keepticks');
datetick( 'y', 'yyyy mmm dd', 'keeplimits');
xtickangle(45);
xlabel('Jupiter departure date [MJD2000]', 'FontSize', 18);
ylabel('Earth flyby date [MJD2000]', 'FontSize', 18);
set(Color,'YtickLabel');
Color.Label.String = '$\Delta V1 [km/s]$';
Color.Label.Interpreter = 'latex';
grid on
title('\textbf{Departure - Flyby Pork Chop Plot}', 'FontSize', 18);

% plot time of flight in the same contour
[C1, h1] = contour(x, y, (T11/(24*3600))','Black');
clabel(C1,h1);                                              

%% PRELIMINARY ANALYSIS OF SECOND LAMBERT ARC
   
arr_window = linspace(departure.date_min, arrival.date_max, temp.n);

temp.Dvmin = astroConstants(5);                             
r = NaN;
c = NaN;                                                  

planes = (1:1:10);                                

for j = 1 : temp.n
    kep2 = uplanet(fb_window(j), 3);
    [r2(:,j), v2(:,j)] = kep2car(kep2(1),kep2(2),kep2(3),kep2(4),kep2(5),kep2(6),astro.muS);

    for k = 1 : temp.n
        kep3 = uplanet(arr_window(k), 2);
        [r3(:,k), v3(:,k)]= kep2car(kep3(1),kep3(2),kep3(3),kep3(4),kep3(5),kep3(6),astro.muS);
        
        T22(j,k)= (arr_window(k)-fb_window(j))*24*60*60;   

        [aT2(j,k),pT,E_T,ERROR,Vi2(j,k,:),Vf2(j,k,:),TPAR,THETA] = ...
            lambertMR(r2(:,j),r3(:,k),T22(j,k),astro.muS,options.orbitType,options.Nrev,options.Ncase,options.LambOptions);
        Dv2x= v3(1,k)-Vf2(j,k,1);
        Dv2y= v3(2,k)-Vf2(j,k,2);
        Dv2z= v3(3,k)-Vf2(j,k,3);
        Dv2= norm([Dv2x Dv2y Dv2z]);   
        
        Dvtot(j,k) = Dv2;    

        if Dvtot(j,k) < temp.Dvmin 
            temp.Dvmin = Dvtot(j,k);
            r = j;
            c = k;
        end
    end
    x(j) = datenum(mjd20002date(fb_window(j)));  
    y(j) = datenum(mjd20002date(arr_window(j))); 
end

% dates obtained from overall problem
fb_day = mjd20002date(results.dates(2));
a_day = mjd20002date(results.dates(3));

pcpL2 = figure();
pcpL2.WindowState = 'maximized';
hold on
[C, h] = contour(x, y, Dvtot', floor(temp.Dvmin) + (planes));
plot(x(r), y(c), '*', 'LineWidth',12);
caxis([temp.Dvmin temp.Dvmin + planes(length(planes))]);                                     
Color = colorbar;    
colormap('Turbo');
clabel(C,h, floor(temp.Dvmin) + (planes));                         
datetick( 'x', 'yyyy mmm dd', 'keeplimits', 'keepticks');
datetick( 'y', 'yyyy mmm dd', 'keeplimits');
xtickangle(45);
xlabel('Earth flyby date [MJD2000]', 'FontSize', 18);
ylabel('Venus arrival date [MJD2000]', 'FontSize', 18);
set(Color,'YtickLabel');
Color.Label.String = '$\Delta V2 [km/s]$';
Color.Label.Interpreter = 'latex';
grid on
title('\textbf{Flyby - Arrival Pork Chop Plot}', 'FontSize', 18);

[C1, h1] = contour(x, y, (T22/(24*3600))','Black');
clabel(C1,h1);

%% SAVING

savefig(pcpL1, '.\Figures\pcpL1.fig');
close (pcpL1)
savefig(pcpL2, '.\Figures\pcpL2.fig');
close (pcpL2)

end

%% OPENING SAVED FIGURES

pcpL1 = openfig('.\Figures\pcpL1.fig');
pcpL1.WindowState = 'maximized';
scatter(datenum(mjd20002date(results.dates(1))), datenum(mjd20002date(results.dates(2))), 'LineWidth', 12);
pcpL2 = openfig('.\Figures\pcpL2.fig');
pcpL2.WindowState = 'maximized';
scatter(datenum(mjd20002date(results.dates(2))), datenum(mjd20002date(results.dates(3))), 'LineWidth', 12);