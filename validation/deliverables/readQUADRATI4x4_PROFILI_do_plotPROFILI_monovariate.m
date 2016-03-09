% readQUADRATI4x4_PROFILI_do_plotPROFILI.m
%
% read the monthly map 1x1 degree x 6 layer
% plot evolution of selected points for 1 year of simulation (put also data

% do plot monovariate
clear all
close all

SI=1;NO=0
curdir=pwd;
Lrun0='10' % OLD RUN 
Lrun1='17'

CFR2RUN=NO % SI =faccio il confronto tra due run; NO non lo faccio
if CFR2RUN==SI
outdir0=['/Users/gcossarini/PROGETTI/COPERNICUS/CALCIFIER/run_0414_0515_carbonatic_' Lrun0 '/'];
end

outdir1=['/Users/gcossarini/PROGETTI/COPERNICUS/CALCIFIER/run_0414_0515_carbonatic_' Lrun1 '/'];
figdir1=['/Users/gcossarini/PROGETTI/COPERNICUS/CALCIFIER/run_0414_0515_carbonatic_' Lrun1 '/figPROFILI/'];



% read AC_ DIC and pHinsitu
VARname='AC_';
nomefile=[outdir1 'PROF_18aree_' VARname '.nc'];
disp(['reading MOD file ' nomefile])
eval(['Mdepth=ncread(''' nomefile ''',''depth'');']); % centro cella
eval(['Marea=ncread(''' nomefile ''',''area'');']);
eval(['Mtempo=ncread(''' nomefile ''',''time'');']);
eval(['MAC_1=ncread(''' nomefile ''',''' VARname ''');']); %[18 x 62 x 15]
MAC_1(MAC_1>1.e19)=NaN;

Mtime=datenum(2014,double(3+Mtempo),15);
VARname='DIC';
nomefile=[outdir1 'PROF_18aree_' VARname '.nc'];
eval(['MDIC1=ncread(''' nomefile ''',''' VARname ''');']); %[18 x 62 x 15]
MDIC1(MDIC1>1.e19)=NaN;

VARname='PH_';
nomefile=[outdir1 'PROF_18aree_' VARname '.nc'];
eval(['MPH_1=ncread(''' nomefile ''',''' VARname ''');']); %[18 x 62 x 15]
MPH_1(MPH_1>1.e19)=NaN;

narea=18;
ndepth=62;

if CFR2RUN==SI
end

% read experimental data in 18 aree of 4 x 4  and depth layers: as in [0 25 50 100:50:200 300:100:1000 1200:200:2500]; 
    LIMDEPTH=1500;
cd('/Users/gcossarini/MEDITERRANEO/dataset_carbsys/')
load QUADRATI4x4_PROFILI_fromDATASET1999_2013.mat 
% it contains: DEPTH [21] CENTRO [20] LAT1 [18] LON1 [18] 
% MAPO3hm MAPO3hs MAPO3cm MAPO3cs MAPPH_insitum MAPPH_insitus MAPALKm MAPALKs MAPDICm ...
% MAPDICs MAPDIC_RICm MAPDIC_RICs MAPPH25m MAPPH25s MAPTm MAPTs MAPSm MAPSs MAPDENSm MAPDENSs 
%  MAPPH_insitu_SWSm MAPPH_insitu_SWSs MAPPH_insitu_0_SWSm MAPPH_insitu_0_SWSs
% [18 x 20]
cd(curdir);
%%
% plot 
%%% PLOT PROFILES
% WHICH VARIABLE HAT TO BE PLOTTED ?????
varname='ALK';PPP=MAC_1;   OBSM=MAPALKm; OBSS=MAPALKs; varxlabel='ALK \mumol/kg';VARlim=[2400 2720];
%varname='DIC';PPP=MDIC1; OBSM=MAPDICm; OBSS=MAPDICs; varxlabel='DIC \mumol/kg';VARlim=[2100 2360];
%varname='PH_';PPP=MPH_1; OBSM=MAPPH_insitu_0_SWSm; OBSS=MAPPH_insitu_0_SWSs; varxlabel='PH@Tinsitu 0dB SWS';VARlim=[7.98 8.2];

for ii=1:1: narea
figure('position',[50 50 400 700]);
ax1=axes('position',[0.07 0.72 0.92, 0.27]); hold on
  plot_MEDcoastline; hold on
  set(gca,'xlim',[-8 36],'ylim',[30 47]);
  for iii=1:1:18
   plot([LON1(iii) LON1(iii)+4 LON1(iii)+4 LON1(iii) LON1(iii)],[LAT1(iii) LAT1(iii) LAT1(iii)+4 LAT1(iii)+4 LAT1(iii)],':','color',[0.5 0.5 0.5],'linewidth',1);
  end
  plot([LON1(ii) LON1(ii)+4 LON1(ii)+4 LON1(ii) LON1(ii)],[LAT1(ii) LAT1(ii) LAT1(ii)+4 LAT1(ii)+4 LAT1(ii)],'-r','linewidth',2);
  set(gca,'fontsize',14);
%  text(-5, 33, ['Lon: ' num2str(LON1(ii)) '-' num2str(LON1(ii)+4) '; Lat: ' num2str(LAT1(ii)) ' - ' num2str(LAT1(ii)+4)],'fontsize',14);         
%  text(-5, 33, ['Lon: ' num2str(LON1(ii)) '-' num2str(LON1(ii)+4) '; Lat: ' num2str(LAT1(ii)) ' - ' num2str(LAT1(ii)+4)],'fontsize',14);         
%  text(-5, 31, ['Run:' Lrun1 '; Apr2014-Jun2015' ],'fontsize',14);         
D4=[1 12 18 26];
C4=repmat([0 0.2 0.4 0.6],3,1);

ax4=axes('position',[0.11 0.07 0.85 0.61]); hold on
 
Mini=4 % start from month Mini
if Mini==1
colori=jet(15);
legendtext={'A','M','J','J','A','S','O','N','D','J05','F','M','A','M','J','OBS'};
elseif Mini==4
colori=jet(12);
legendtext={'J';'A';'S';'O';'N';'D';'J05';'F';'M';'A';'M';'J';'OBS'};

end

       for im=Mini:1:15
           [p(im-Mini+1)]=plot(PPP(ii,:,im),-Mdepth,'-','color',colori(im-Mini+1,:),'linewidth',2); hold on
       end
       [p2]=plot(squeeze(OBSM(ii,:)), -CENTRO, 'o-k','markersize',5,'linewidth',1.5); hold on
       plot(squeeze(OBSM(ii,:)-OBSS(ii,:)), -CENTRO, '--k','linewidth',1.5);
       plot(squeeze(OBSM(ii,:)+OBSS(ii,:)), -CENTRO, '--k','linewidth',1.5);
       [x1]=xlabel(varxlabel);set(x1,'fontsize',16); set(gca,'fontsize',14);
       set(gca,'ylim',[-LIMDEPTH 0],'xlim',VARlim); grid on
 %      L1=legend([p p2],'A|M|J|J|A|S|O|N|D|J05|F|M|A|M|J|OBS',0);set(L1,'fontsize',14);
       L1=legend([p p2],legendtext,0);set(L1,'fontsize',14);
box on
  
set(gcf,'paperpositionmode','auto','renderer','zbuffer')

if ii<10
nomefileout=[figdir1 'Q4x4_PROF_' varname 'runs'  Lrun1 '_A0' num2str(ii) '.png'];
else
nomefileout=[figdir1 'Q4x4_PROF_' varname 'runs'  Lrun1 '_A' num2str(ii) '.png'];
end
eval (['print ' nomefileout ' -dpng -r300' ]);

% calcolo le statistiche di errore
OO=squeeze(OBSM(ii,:));
MM=squeeze(mean(PPP(ii,:,Mini:end),3));
m=interp1(Mdepth,MM,CENTRO,'linear');
o=OO(:);
m=m(:);

indnan=(isnan(o) | isnan(m) );
o(indnan)=[];
m(indnan)=[];

d=m-o;
bias(ii,1)=mean(d);
rmsd(ii,1)=(mean(d.^2))^.5;
c=corrcoef(m,o); corr(ii,1)=c(2,1);
end
%%
bias=double(bias);
corr=double(corr);
rmsd=double(rmsd);
nomefile= [figdir1 'RUN' Lrun1 '_' varname '_bias.txt'];
bias(isnan(bias))=-999.;
save(nomefile,'bias','-ascii','-tabs')
nomefile= [figdir1 'RUN' Lrun1 '_' varname 'rmsd.txt'];
rmsd(isnan(rmsd))=-999.;
 save(nomefile,'rmsd','-ascii','-tabs')
nomefile= [figdir1 'RUN' Lrun1 '_' varname '_corr.txt'];
corr(isnan(corr))=-999.;
 save(nomefile,'corr','-ascii','-tabs')



 

%%
% read dataset carbosys
cd ('/Users/gcossarini/MEDITERRANEO/dataset_carbsys/')
load DATASET_MED_CARBSYS1999_2013.mat
cd(curdir)

ind=find(DATI(:,6)<20);
 % plot di tutti i quadrati
 figure('position',[50 50 600 400]);
ax1=axes('position',[0.1 0.1 0.89, 0.89]); hold on
  plot_MEDcoastline; hold on
  set(gca,'xlim',[-8 36],'ylim',[30 47]);
  plot(DATI(ind,5),DATI(ind,4),'.','markersize',8);
 for iii=1:1:18
   plot([LON1(iii) LON1(iii)+4 LON1(iii)+4 LON1(iii) LON1(iii)],[LAT1(iii) LAT1(iii) LAT1(iii)+4 LAT1(iii)+4 LAT1(iii)],'-','color',[0.5 0.5 0.5],'linewidth',2);
   text(LON1(iii)+2,LAT1(iii)+2,num2str(iii),'fontsize',20,'verticalalignment','middle','horizontalalignment','center');
  end
  box on
  set(gca,'fontsize',14);
set(gcf,'paperpositionmode','auto','renderer','zbuffer')
nomefile= [figdir1 'MAP_MED_quadrati4x4.png'];
eval (['print ' nomefile ' -dpng -r300' ]);

