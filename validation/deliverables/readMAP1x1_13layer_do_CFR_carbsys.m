% readMAP1x1_6layer_do_CFR_carbsys.m
%
% read the monthly map 1x1 degree x 6 layer
% compare with the data set carbsys2

clear all
%close all

curdir=pwd;

NLAYER=7
outdir='/Users/gcossarini/PROGETTI/COPERNICUS/CALCIFIER/run_0414_0515/'
figdir='/Users/gcossarini/PROGETTI/COPERNICUS/CALCIFIER/run_0414_0515/figMAPPE/'

runname='17';
outdir='/Users/gcossarini/PROGETTI/COPERNICUS/CALCIFIER/run_0414_0515_carbonatic_02/'
figdir='/Users/gcossarini/PROGETTI/COPERNICUS/CALCIFIER/run_0414_0515_carbonatic_02/figMAPPE/'
outdir='/Users/gcossarini/PROGETTI/COPERNICUS/CALCIFIER/run_0414_0515_carbonatic_03/'
figdir='/Users/gcossarini/PROGETTI/COPERNICUS/CALCIFIER/run_0414_0515_carbonatic_03/figMAPPE/'

outdir='/Users/gcossarini/PROGETTI/COPERNICUS/CALCIFIER/run_0414_0515_carbonatic_05/'
figdir='/Users/gcossarini/PROGETTI/COPERNICUS/CALCIFIER/run_0414_0515_carbonatic_05/figMAPPE/'
outdir=['/Users/gcossarini/PROGETTI/COPERNICUS/CALCIFIER/run_0414_0515_carbonatic_' runname '/']
figdir=['/Users/gcossarini/PROGETTI/COPERNICUS/CALCIFIER/run_0414_0515_carbonatic_' runname '/figMAPPE/']

VARname='DIC'; VARUNIT='\mumol/kg'; VARLIM=[2050 2400]; VARdati='DIC'; ERRLIM=[-80 80];
%VARname='O3c'; VARUNIT='mgC/m^3'; VARLIM=[26000 29000]; VARdati='O3C'; ERRLIM=[-1000 1000];
%VARname='O3h'; VARUNIT='mmol/m^3'; VARLIM=[2400 2800]; VARdati='O3H'; ERRLIM=[-80 80];
VARname='AC_'; VARUNIT='\mumol/kg'; VARLIM=[2380 2750]; VARdati='ALK'; ERRLIM=[-80 80];
%VARname='PH_'; VARUNIT=' '; VARLIM=[8.0 8.2]; VARdati='PH_SWS_Tins_0dB'; ERRLIM=[-0.05 0.05];
%VARname='pCO'; VARUNIT=' '; VARLIM=[300 420]; VARdati='PCO'; ERRLIM=[-0.05 0.05];
%VARname='ppn'; VARUNIT='mgC/m3/d'; VARLIM=[-1 10]; VARdati=''; ERRLIM=[-5 5];
%VARname='N1p'; VARUNIT='\mumol/m^3'; VARLIM=[0 0.25]; VARdati=''; ERRLIM=[-0.1 0.1];
%VARname='N3n'; VARUNIT='\mumol/m^3'; VARLIM=[0 6]; VARdati=''; ERRLIM=[-0.1 0.1];
%VARname='R6c'; VARUNIT='mgC/m^3'; VARLIM=[0 12]; VARdati=''; ERRLIM=[-0.1 0.1];
%VARname='R2c'; VARUNIT='mgC/m^3'; VARLIM=[0 30]; VARdati=''; ERRLIM=[-0.1 0.1];


% dimension of Model file 43 x 16 x 13 x 15
Depth1=[00, 050, 100, 200, 0500, 1500];
Depth2=[50, 100, 200, 500, 1500, 4000];
STRATIname={'0-50','50-100','100-150', '150-200','200-500','500-1000', '1000-1500','1500-2000', '2000-2500', '2500-3000','3000-3500','3500-4000','4000-4500'};
nlat=17;
nlon=44;
nlayer=13;
ntime=15;
nomefile=[outdir 'MAP1x1_13lev_15m_' VARname '.nc'];
disp(['reading MOD file ' nomefile])
% file MONTHLY MAPs [15 frames] from 04/2014 to 06/2015
eval(['Mlayer=ncread(''' nomefile ''',''depth'');']); % centro cella
eval(['MLon=ncread(''' nomefile ''',''lon'');']);
eval(['MLat=ncread(''' nomefile ''',''lat'');']);
eval(['MOD=ncread(''' nomefile ''',''' VARname ''');']); %[43 x 16 x 6 x 15]
[MXLON, MYLAT]=meshgrid([MLon; MLon(end)+1],[MLat; MLat(end)+1]); % ripeto ultimo punto per il plot

MOD(MOD>1E19)=NaN;
MXLON=double(MXLON); MYLAT=double(MYLAT); MOD=double(MOD);

    

% read experimental data in 1 x 1  and 6 depth
OBSlayername={'0-50','50-100','100-150','150-200','200-500','500-1500','1500-4000'}; 
OBSlev=[25,75,125,175,350,1000,2500];
cd('/Users/gcossarini/MEDITERRANEO/dataset_carbsys/')
disp(['reading OBS file '  VARdati '_lon_lat_1GRADO_7depth_5stag.mat'])
if exist([VARdati '_lon_lat_1GRADO_7depth_5stag.mat' ], 'file') == 2
    eval(['load  ' VARdati '_lon_lat_1GRADO_7depth_5stag.mat '])
     eval(['OBS=LONLATDEPTHSTAG_' VARdati ';']);
else
    disp('file not existing')
    OBS(1:44,1:17,1:7,1:5,1:3)=NaN;
    Xlon=[-8:1:35];Ylat=[30:1:46];
end
% LONLATDEPTHSTAG_ALK Xlon Ylat STRATIname nomestag -MAT']);
cd(curdir); %  [44 17 7(strati) 5(stag) 3(stat: mean, std, #, median, p25, p75)]

[OXLON, OYLAT]=meshgrid([Xlon, Xlon(end)+1], [Ylat, Ylat(end)+1]); % aggiungo un punto per il pcolor

% MOD da 04/2014 a 06
% MEDIA annaule del modello DA 07-2014 a 06-2015 --> da 4 a 15

% CHOICE WHICH MODEL FRAME TO USE
MODm=mean(MOD(:,:,:,4:15),4); labeltime=' ave 07/14 - 06/15'; labelfig='aveYEAR'; % annual mean
%MODm=mean(MOD(:,:,:,3),4); labeltime='  06/14 '; labelfig='201406'; % annual mean
%MODm=mean(MOD(:,:,:,15),4); labeltime='  06/15 '; labelfig='201506'; % annual mean

NLEV=7;
%NLEV=1;

Pylim=[30 46.5];
Pxlim=[-6 36];
for ilev=1:1:NLEV
    disp(['plotting level ' num2str(ilev)]);
figure('position',[50 50 600 800]); set(gcf,'renderer','zbuffer','paperpositionmode','auto');
subplot(3,1,1);
if ilev<6
MM=squeeze(MODm(:,:,ilev)); 
elseif ilev==6
MM=squeeze(nanmean(MODm(:,:,6:7),3)); % 500-1000 + 1000-1500 
elseif ilev==7
MM=squeeze(nanmean(MODm(:,:,8:12),3)); % 1500-2000 + 2000-2500 +  2500-3000 + 3000-3500 +  3500-4000 
end    
PMM=MM; PMM(:,end+1)=PMM(:,end);PMM(end+1,:)=PMM(end,:);

pcolor(MXLON-.5,MYLAT-.5,PMM'); shading flat; hold on; caxis(VARLIM); [c1]=colorbar; set(c1,'fontsize',12);
plot_MEDcoastline
set(gca,'xlim',Pxlim,'ylim',Pylim,'fontsize',12);
T1=title([VARname ' ::: MODEL ' labeltime ' -  lev: ' OBSlayername{ilev}] );set(T1,'fontsize',15);

subplot(3,1,2);

OO=squeeze(OBS(:,:,ilev,1));
% verifico che i dati non siano troppo pochi
%% TBD
%


POO=OO; POO(:,end+1)=POO(:,end);POO(end+1,:)=POO(end,:);
pcolor(OXLON-.5,OYLAT-.5,POO'); shading flat; hold on; caxis(VARLIM); [c1]=colorbar; set(c1,'fontsize',12);
plot_MEDcoastline
set(gca,'xlim',Pxlim,'ylim',Pylim,'fontsize',12);
T1=title('OBS');set(T1,'fontsize',14);
subplot(3,1,3)
ERR=MM-OO;
PERR=ERR; PERR(:,end+1)=PERR(:,end);PERR(end+1,:)=PERR(end,:);

pcolor(OXLON-.5,OYLAT-.5,PERR'); shading flat; hold on; caxis(ERRLIM); [c1]=colorbar; set(c1,'fontsize',12);
plot_MEDcoastline
set(gca,'xlim',Pxlim,'ylim',Pylim,'fontsize',12);
T1=title('MOD-OBS');set(T1,'fontsize',14);

nomefileout=[figdir 'CFR_M_O_MAP1x1_' VARname '_lev' num2str(ilev) '_'  labelfig '.png'];
%eval (['print ' nomefileout ' -dpng -r300' ]);

%% calcolo le statistiche di errore

for ibac=1:1:1
    OO(1:3,:)=NaN; % not use data on the cadiz bay (atlantic buffer)
o=OO(:);
m=MM(:);

indnan=(isnan(o) | isnan(m) );
o(indnan)=[];
m(indnan)=[];

d=m-o;
if VARname=='DIC';
     d(d>50)=[]; d(d<-50)=[]; % elimino gli errori + alti
end
if VARname=='AC_';
     d(d>60)=[]; d(d<-60)=[]; % elimino gli errori + alti
end

bias(ilev,ibac)=mean(d);
rmsd(ilev,ibac)=(mean(d.^2))^.5;
c=corrcoef(m,o); corr(ilev,ibac)=c(2,1);
end
end

break

nomefile= [figdir 'RUN' runname '_' VARname '_bias.txt'];
 save(nomefile,'bias','-ascii','-tabs')
nomefile= [figdir 'RUN' runname '_' VARname 'rmsd.txt'];
 save(nomefile,'rmsd','-ascii','-tabs')
nomefile= [figdir 'RUN' runname '_' VARname '_corr.txt'];
 save(nomefile,'corr','-ascii','-tabs')


