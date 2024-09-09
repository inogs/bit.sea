for idate in `ls STATS/*orig*ORIG* | xargs -n 1 basename| cut -c 11-18 `
 do
  echo mv STATS/stats_orig${idate}ORIG.pkl STATS_time${idate}ORIG.pkl
done
