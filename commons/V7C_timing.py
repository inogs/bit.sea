from __future__ import print_function
from datetime import datetime,timedelta


def last_day(d, day_name):
    '''
    Takes in account that at 09:00 V6C has analysis in archive '''
    days_of_week = ['sunday','monday','tuesday','wednesday',
                        'thursday','friday','saturday']
    target_day = days_of_week.index(day_name.lower())
    delta_day = target_day - d.isoweekday()
    if delta_day > 0: delta_day -= 7 # go back 7 days
    if day_name.lower()=='tuesday':
        if delta_day ==0 :
            if d.hour < 9:
                delta_day -= 7
            else:
                target_day - d.isoweekday()
    return d + timedelta(days=delta_day)

def next_day(d, day_name):
    days_of_week = ['sunday','monday','tuesday','wednesday',
                        'thursday','friday','saturday']
    target_day = days_of_week.index(day_name.lower())
    delta_day = target_day - d.isoweekday()
    if delta_day <= 0: delta_day += 7 # go back 7 days
    return d + timedelta(days=delta_day)

def round(timeobj):
    return datetime.strptime(timeobj.strftime("%Y%m%d"), "%Y%m%d")

today   =  datetime.now()
yesterday=  round(today - timedelta(days=1))
Last_rundate_analysis=round(last_day(today, 'tuesday'))
date__end=Last_rundate_analysis - timedelta(days=1)


def find_best(d):
    '''
     V6C archive directory for a specific date
    The corresponding ave files in that directory are the best choice for that date,
    and correspond to DU content.
    Arguments:
    * d * datetime object
    Returns:
    * mtype   * string 'Analysis', 'Hindcast' or 'Forecast'
    * rundate * a datetime object 
    '''    
    if d<=date__end: #analysis
        if d.isoweekday() == 1:
            Analysis_rundate=next_day(d+timedelta(days=7), "tuesday")
        else:
            Analysis_rundate=next_day(d, "tuesday")
        mtype="Analysis"
        rundate=Analysis_rundate

    else:
        # searching in forecast
        if d<=yesterday:
            mtype="Hindcast"
            rundate=d+timedelta(days=1)
        else:
            mtype="Forecast"
            rundate=yesterday
    return mtype, rundate

def find_best_dir(datestr, dateformat="%Y%m%d"):
    '''Prints V6C archive directory for a specific date
    The corresponding ave files in that directory are the best choice for that date,
    and correspond to DU content.
    
    It prints something like that
    analysis/20200505
    forecast/20200507'''
        
    d=datetime.strptime(datestr,dateformat)
    if d.hour==0: d += timedelta(hours=12)
    mtype, rundate= find_best(d)
    if mtype in ['Hindcast',"Forecast"]:
        runtype="forecast"
    else:
        runtype="analysis"
    return runtype + "/" + rundate.strftime("%Y%m%d")

def letter(mtype):
    if mtype=='Hindcast' : return 's'
    if mtype=='Forecast' : return 'f'
    if mtype=='Analysis' : return 'a'

def find_best_forcing(datestr,dateformat="%Y%m%d"):
    d=datetime.strptime(datestr,dateformat)
    if d.hour==0: d += timedelta(hours=12)
    mtype, rundate= find_best(d)
    thedir=find_best_dir(datestr)
    forcing="%s/CMCC_PHYS/mfs_eas6-%s-%s-%s-"  %(thedir, rundate.strftime("%Y%m%d"), d.strftime("%Y%m%d"), letter(mtype))
    return forcing
def find_best_bgc(datestr,dateformat="%Y%m%d"):
    thedir=find_best_dir(datestr)
    return thedir + "/POSTPROC/AVE_FREQ_1/ARCHIVE/ave." + datestr + "-12:00:00."

def list_for_maps(datestr,dateformat="%Y%m%d"):
    timeobj=datetime.strptime(datestr,dateformat)
    DIRS=[]
    for deltadays in range(-10,11):
        d=timeobj + timedelta(days=deltadays)
        thedir=find_best_dir(d.strftime("%Y%m%d"))
        if thedir not in DIRS:
            DIRS.append(thedir)
    return DIRS



if __name__=="__main__":

    print(list_for_maps("20200512"))

    from commons import genUserDateList as DL
    DAYS=DL.getTimeList("20200505-12:00:00", "20200525-12:00:00", "days=1")

    for d in DAYS:
        weekday=d.isoweekday()
        mtype, rundate= find_best(d)
        datestr=d.strftime("%Y%m%d")

        a= "%s %s %s from run of %s"  %(datestr , weekday, mtype,  rundate.strftime("%Y%m%d"))
        forcingT=find_best_forcing(datestr)
        print(a + " " + forcingT)
        

