from datetime import datetime


class season:
    """
    Season Object

    This object provide a tools for season definition and management
    The default setting is classical 4 season, the star day is set as:
    - winter : 1221
    - spring : 0321
    - summer : 0621
    - fall   : 0921
    """


    def __init__(self):
        """
        Set classical for season and init object attribute
        """
        self.numbers_season = 0
        self.SEASON_LIST = []
        self.setseasons(["1221","0321","0622","0921"])

    def setseasons(self,stendseason):
        """
        Given an arry where is defined the date when season start and generate
        the season list. In input take an array of string.
        """
        self.SEASON_LIST = []
        self.numbers_season = len(stendseason)
        for i in stendseason:
            self.SEASON_LIST.append(datetime.strptime(i,'%m%d'))

    def get_seasons_number(self):
        """
        Return the number of seasons defined in object
        """
        return self.numbers_season

    def get_season_dates(self,season_num):
        """
        Given season number, return the date of start and end.
        """
        season = []
        if (season_num >= self.numbers_season):
            print("Numbers of seasons is less than",season_num)
        season.append(self.SEASON_LIST[season_num])

        season_num = season_num + 1
        if (season_num >= self.numbers_season):
            season_num = 0

        season.append(self.SEASON_LIST[season_num])

        return season

    def findseason(self,date):
        """
        Take a date as input and return the season where it is in
        """

        check = -1
        day   = ( date.month, date.day )
        for i in range(0,self.numbers_season):
            season = self.get_season_dates(i)
            #print season[0],season[1],"\n"
            #print date
            begin = ( season[0].month, season[0].day )
            end   = ( season[1].month, season[1].day )

            if( begin <= day < end):
                check = i

            #print check

        if(check == -1):
            season = self.get_season_dates(0)
            begin   = ( season[0].month, season[0].day )
            end   = ( season[1].month, season[1].day )
            if(  day == begin):
                check = 0
            if(  day < end):
                check = 0
            if(  day == end):
                check = 1

        if(check == -1):
            season = self.get_season_dates(self.numbers_season-1)
            begin = ( season[0].month, season[0].day )
            end   = ( season[1].month, season[1].day )
            if(  day > end):
                check = 0

        assert (check >= 0),"ERROR : date case in not contepleted! :-( "

        return check
