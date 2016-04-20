from datetime import datetime

#author : epascolo

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
        Init object and set astronomical season
        """
        self.numbers_season = 0
        self.SEASON_LIST = []
        self.SEASON_LIST_NAME = []
        self.setseasons(["1221","0321","0622","0921"],["winter","spring","summer","fall"])

    def setseasons(self,startseason,nameseason):
        """

        Given two arrays where is defined the date when season start and its name
        the subroutine generate the season list. In input take two arrays of string:

        - Example of startseason : ["1221","0321","0622","0921"]
        - Example of nameseason  : ["winter","spring","summer","fall"]

        In the previous example we shown that winter start to 21 december,
        spring to 21 march, summer to 22 june and fall on 21 september.

        """
        self.SEASON_LIST = []
        self.SEASON_LIST_NAME = []
        self.numbers_season = len(startseason)

        if (len(startseason) != len(nameseason)):
            print("ERROR : arrays definitions mistmatch! :-( ")
            exit()

        for i in range(0,self.numbers_season):
            self.SEASON_LIST.append(datetime.strptime(startseason[i],'%m%d'))
            self.SEASON_LIST_NAME.append(nameseason[i])

    def get_seasons_number(self):
        """
        Return the number of seasons defined in this object
        """
        return self.numbers_season

    def get_season_dates(self,season_num):
        """
        Given season number, return the range of season dates (start and end)
        and the name of season.
        """
        season = []
        if (season_num >= self.numbers_season):
            print("Numbers of seasons is less than",season_num)
        season.append(self.SEASON_LIST[season_num])

        season_num = season_num + 1
        if (season_num >= self.numbers_season):
            season_num = 0

        season.append(self.SEASON_LIST[season_num])
        season_name = self.SEASON_LIST_NAME[season_num]

        return season,season_name

    def findseason(self,date):
        """
        Takes a date as input and return the number and name of season where it is in.

        """

        check = -1
        check_name = "Not reconized"
        day   = ( date.month, date.day )
        for i in range(0,self.numbers_season):
            season,name = self.get_season_dates(i)
            #print season[0],season[1],"\n"
            #print date
            begin = ( season[0].month, season[0].day )
            end   = ( season[1].month, season[1].day )

            if( begin <= day < end):
                check = i
                check_name = name

            #print check

        if(check == -1):
            season,name = self.get_season_dates(0)
            begin   = ( season[0].month, season[0].day )
            end   = ( season[1].month, season[1].day )
            if(  day == begin):
                check = 0
                check_name = self.SEASON_LIST_NAME[0]
            if(  day < end):
                check = 0
                check_name = self.SEASON_LIST_NAME[0]
            if(  day == end):
                check = 1
                check_name = self.SEASON_LIST_NAME[1]

        if(check == -1):
            season,name = self.get_season_dates(self.numbers_season-1)
            begin = ( season[0].month, season[0].day )
            end   = ( season[1].month, season[1].day )
            if(  day > end):
                check = 0
                check_name = name

        assert (check >= 0),"ERROR : date case in not manged! :-( "

        return check,check_name
