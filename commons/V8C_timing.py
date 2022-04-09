from __future__ import print_function
from datetime import datetime,timedelta

def timeround(timeobj):
    return datetime.strptime(timeobj.strftime("%Y%m%d"), "%Y%m%d")

def last_analysis_rundate(d, analysis_day='tuesday', sure_an_already_run=False):
    '''
    Takes in account that at 09:00 V6C has analysis in archive '''
    days_of_week = ['sunday','monday','tuesday','wednesday',
                        'thursday','friday','saturday']
    target_day = days_of_week.index(analysis_day.lower())
    delta_day = target_day - d.isoweekday()
    if delta_day > 0: delta_day -= 7 # go back 7 days
    if analysis_day.lower()=='tuesday':
        if delta_day ==0 :
            if sure_an_already_run:
                delta_day = target_day - d.isoweekday()
            else:
                if d.hour < 9:
                    delta_day -= 7
                else:
                    delta_day = target_day - d.isoweekday()
    return timeround(d + timedelta(days=delta_day))

def next_day(d, day_name):
    days_of_week = ['sunday','monday','tuesday','wednesday',
                        'thursday','friday','saturday']
    target_day = days_of_week.index(day_name.lower())
    delta_day = target_day - d.isoweekday()
    if delta_day <= 0: delta_day += 7 # go back 7 days
    return d + timedelta(days=delta_day)


def get_last_analysis_day(sure_an_already_run=False, today=datetime.now()):
    '''
    Argument:
    * sure_an_already_run * logical, true if we are sure that last analysis is already archived
    Returns
    * last_analysis_day * datetime object
    '''
    Last_rundate_analysis=last_analysis_rundate(today, sure_an_already_run=sure_an_already_run)
    return Last_rundate_analysis - timedelta(days=1)

def get_last_forecast_rundate(sure_fc_already_run=False, today=datetime.now()):
    '''
    Argument:
    * sure_fc_already_run * logical, true if we are sure that last forecast is already archived
    Returns
    * last_forecast rundate * datetime object
    '''
    if sure_fc_already_run:
        last_forecast_rundate=timeround(today)
    else:
        last_forecast_rundate =  timeround(today - timedelta(days=1))
    return last_forecast_rundate



def find_best(d, analisis_end_date, last_forecast_rundate):
    '''
     V8C archive directory for a specific date
    The corresponding ave files in that directory are the best choice for that date,
    and correspond to DU content.
    Arguments:
    * d * datetime object
    Returns:
    * mtype   * string 'Analysis', 'Hindcast' or 'Forecast'
    * rundate * a datetime object 
    '''    
    if d<=analisis_end_date: #analysis
        if d.isoweekday() == 1:
            Analysis_rundate=next_day(d+timedelta(days=7), "tuesday")
        else:
            Analysis_rundate=next_day(d, "tuesday")
        mtype="Analysis"
        rundate=Analysis_rundate

    else:
        # searching in forecast
        if d<=last_forecast_rundate:
            mtype="Hindcast"
            rundate=d+timedelta(days=1)
        else:
            mtype="Forecast"
            rundate=last_forecast_rundate
    return mtype, rundate

def find_best_dir(datestr, dateformat="%Y%m%d", sure_an_already_run=True, sure_fc_already_run=True):
    '''Prints V6C archive directory for a specific date
    The corresponding ave files in that directory are the best choice for that date,
    and correspond to DU content.
    
    It prints something like that
    analysis/20200505
    forecast/20200507'''
    
    analysis_end_date = get_last_analysis_day(sure_an_already_run)
    last_forecast_rundate=get_last_forecast_rundate(sure_fc_already_run)

    d=datetime.strptime(datestr,dateformat)
    if d.hour==0: d += timedelta(hours=12)
    mtype, rundate= find_best(d,analysis_end_date,last_forecast_rundate)
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
    mtype, rundate= find_best(d,analysis_end_date)
    thedir=find_best_dir(datestr)
    forcing="%s/CMCC_PHYS/mfs_eas6v8-%s-%s-%s-"  %(thedir, rundate.strftime("%Y%m%d"), d.strftime("%Y%m%d"), letter(mtype))
    return forcing
def find_best_bgc(datestr,dateformat="%Y%m%d"):
    thedir=find_best_dir(datestr)
    return thedir + "/POSTPROC/AVE_FREQ_1/ARCHIVE/ave." + datestr + "-12:00:00."

def list_for_maps(datestr,dateformat="%Y%m%d",sure_an_already_run=False, sure_fc_already_run=False):
    timeobj=datetime.strptime(datestr,dateformat)
    DIRS=[]
    for deltadays in range(-10,11):
        d=timeobj + timedelta(days=deltadays)
        thedir=find_best_dir(d.strftime("%Y%m%d"), sure_an_already_run=sure_an_already_run, sure_fc_already_run=sure_fc_already_run)
        #print(d, thedir)
        if thedir not in DIRS:
            DIRS.append(thedir)
    return DIRS



if __name__=="__main__":
   
    d = datetime(2022,4,5)
    analysis_end_date = get_last_analysis_day(sure_an_already_run=True, today=datetime(2022,4,6,5,12))
    last_fc_date = get_last_analysis_day()

    
    print(analysis_end_date)
    print (find_best(d,analysis_end_date,last_fc_date))

    print(list_for_maps("20220406", sure_an_already_run=False, sure_fc_already_run=True))

    import sys
    sys.exit()

    from commons import genUserDateList as DL
    DAYS=DL.getTimeList("20200505-12:00:00", "20200525-12:00:00", "days=1")

    for d in DAYS:
        weekday=d.isoweekday()
        mtype, rundate= find_best(d,analysis_end_date)
        datestr=d.strftime("%Y%m%d")

        a= "%s %s %s from run of %s"  %(datestr , weekday, mtype,  rundate.strftime("%Y%m%d"))
        forcingT=find_best_forcing(datestr)
        print(a + " " + forcingT)
        

