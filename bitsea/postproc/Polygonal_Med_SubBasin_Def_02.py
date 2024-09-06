class Subbasin():
    def __init__(self,isdef=False,x=0,y=0):
        self.x=x
        self.y=y
        self.isdef = isdef

class Composed_Subbasin():
        def __init__(self,LIST):
            self.isdef=False
            self.LIST=LIST

alb=Subbasin()
alb.x=[-5.5,-1.0,-1.0,-5.5];
alb.y=[32.0,32.0,40.0,40.0];
alb.isdef = True; 

sww=Subbasin()
sww.x=[-1.0, 5,  5,-1.0 ];
sww.y=[32.0,32.00, 39.50,39.5 ];
sww.isdef = True; 

swe=Subbasin()
swe.x=[5, 9.25,  9.25, 5 ];
swe.y=[32.0,32.00, 39.50,39.5 ];
swe.isdef = True; 

nwm=Subbasin()
nwm.x=[-1.0,9.25, 9.25,15.00,13.0, 10,-1];
nwm.y=[39.5,39.5,41.25,41.25,42.5 ,46,46];
nwm.isdef = True; 

tyr=Subbasin()
tyr.x=[ 9.25,15.00,16.5,16.5,15.85,15.85,15.46875,14.6,12.5,12.5, 9.25];
tyr.y=[41.25,41.25,40.0,38.6,38.20,38.30,38.21875,38.0,37.8,36.8,36.80] ;
tyr.isdef = True ; 

adn=Subbasin()
adn.x=[10,13.0,20.0,15];
adn.y=[46,42.5,42.5,46];
adn.isdef = True ; 

ads=Subbasin()
ads.x=[ 14.0,20.0,20,18.5,18.0,16.6,13.0];
ads.y=[ 42.5,42.5,40,40.0,40.5,41.0,42.5];
ads.isdef = True ;

aeg=Subbasin()
aeg.x=[22.80,23,23,   28,   28,22,22.8 ];
aeg.y=[37.15,37,35.25,35.25,43,43,37.15];
aeg.isdef = True ;

ion=Subbasin()
ion.x=[16.5,16.5,15.85,15.85,15.46875,14.6,12.5,12.5,9.25,9.25,21.63,21.63,23.0,20,18.5,18.0,16.6];
ion.y=[40.0,38.6,38.20,38.30,38.21875,38.0,37.8,36.8,36.80,30.0,30.00,36.82,38.7,40,40.0,40.5,41.0];
ion.isdef = True ;

lev=Subbasin()
lev.x=[21.63,21.63,21.90,22.80,23,23,   28 ,  28,36.3,36.3]; 
lev.y=[30.00,36.82,37.15,37.15,37,35.25,35.25,38,38,  30  ];
lev.isdef = True ; 

atl=Subbasin()
atl.x =[ -5.5,-5.5,-9.0,-9.0]; 
atl.y =[ 32.0,40.0,40.0,32.0]; 
atl.isdef = True; 


wes=Composed_Subbasin(['alb', 'sww', 'swe' , 'nwm', 'tyr' ])
eas=Composed_Subbasin(['ion', 'lev' ])
med=Composed_Subbasin(['alb', 'sww', 'swe' , 'nwm', 'tyr', 'ion', 'lev' ])

P=[ 'alb' ,'sww','swe','nwm','tyr','adn','ads','aeg','ion','lev']                            