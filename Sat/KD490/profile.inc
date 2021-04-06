function my_prex_or_die {
  # Like my_prex, but in case of failure runs my_die
  typeset _command="$1"
  shift 1
  my_prex "$_command" || {
    typeset -i _ec=${1:-1}
    typeset    _message="${2:-Command '$_command' failed}"
    my_die "$_ec" "$_message"
  }
}


function my_date {
 date
}

function my_prex {
  # PRint a command and EXecutes it;
  # after the command has been executed, my_prex prints its exit_code and returns the same value
  typeset -i _ec
  printf "PrEx $(my_date)    > %s\n" "$*"
  _logical_sim="${2:-false}"

  if $_logical_sim ; then
    echo "[WRN: CMD SIMULATED]:: $1" ; _ec=$?
  else
    eval "$1" ; _ec=$?
  fi
  typeset _status
  if [[ $_ec -eq 0 ]] ; then
    _status="ok"
  else
    _status="KO"
  fi
  printf "PrEx $(my_date) %3d< %s [%s]\n" "$_ec" "$*" "$_status"
  return $_ec
}

function my_die {
  # log an error message and die
  typeset -i _ec="$1"
  shift 1
  my_log 0 "$@"
  my_exit $_ec
}

function my_log {
  # Write a log messages ($2...) on standard output;
  # is $1 is lower than $my_LOGMAX, the messages is written also onto $my_LOGFILE
  if [[ $my_LOGMAX -eq $my_DISABLE_LOG ]] ; then
    return 0
  fi
  typeset -i _nlog="$1"
  shift 1
  for _message in "$@" ; do
    if ( $my_BATCH ) && [[ $my_LOGMAX -ne $my_DISABLE_LOG_ON_LOGFILE && $_nlog -le $my_LOGMAX ]] ; then
      printf '%s %02d (%s) : %s\n' "$(my_date)" "$_nlog" "$my_RUNID" "$_message" | tee -a $my_LOGFILE
    else
      printf '%s %02d (%s) : %s\n' "$(my_date)" "$_nlog" "$my_RUNID" "$_message"
    fi
  done
}

function my_exit {
  # Log message for script exit
  typeset -i _default_ec=${_my_EXIT_CODE:-$?}
  typeset -i _ec=${1:-$_default_ec}
  #echo "my_exit(${1}->${_default_ec}->${_ec})"
  typeset    _message=""
  if [[ "$1" == "" ]] ; then
    _message=" (from last executed command)"
  fi
  if $_my_START_DONE && ! $_my_EXIT_DONE ; then
    my_log 0 "* EXIT(${my_PROGNAME}) [$(my_exitcode ${_ec})]$_message"
    _my_EXIT_DONE=true
    _my_EXIT_CODE=${_ec}
  fi
  exit $_ec
}


