from datetime import datetime,timedelta

def last_day(d, day_name):
    days_of_week = ['sunday','monday','tuesday','wednesday',
                        'thursday','friday','saturday']
    target_day = days_of_week.index(day_name.lower())
    delta_day = target_day - d.isoweekday()
    if delta_day >= 0: delta_day -= 7 # go back 7 days
    return d + timedelta(days=delta_day)

def next_day(d, day_name):
    days_of_week = ['sunday','monday','tuesday','wednesday',
                        'thursday','friday','saturday']
    target_day = days_of_week.index(day_name.lower())
    delta_day = target_day - d.isoweekday()
    if delta_day <= 0: delta_day += 7 # go back 7 days
    return d + timedelta(days=delta_day)



today   =  datetime.strptime(datetime.now().strftime("%Y%m%d"), "%Y%m%d")
yesterday=  today - timedelta(days=1)
Last_rundate_analysis=last_day(today, 'tuesday')
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
            rundate=d
        else:
            mtype="Forecast"
            rundate=yesterday
    return mtype, rundate

def find_best_datestr(datestr, dateformat="%Y%m%d"):
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
                


if __name__=="__main__":
    from commons import genUserDateList as DL
    DAYS=DL.getTimeList("20200401-12:00:00", "20200510-12:00:00", "days=1")

    for d in DAYS:
        weekday=d.isoweekday()
        mtype, rundate= find_best(d)
        datestr=d.strftime("%Y%m%d")
        print "%s %s %s from run of %s"  %(datestr , weekday, mtype,  rundate.strftime("%Y%m%d"))

        

