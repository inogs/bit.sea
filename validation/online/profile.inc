function opa_prex_or_die {
  # Like opa_prex, but in case of failure runs opa_die
  typeset _command="$1"
  shift 1
  opa_prex "$_command" || {
    typeset -i _ec=${1:-1}
    typeset    _message="${2:-Command '$_command' failed}"
    opa_die "$_ec" "$_message"
  }
}


function opa_date {
 date
}

function opa_prex {
  # PRint a command and EXecutes it;
  # after the command has been executed, opa_prex prints its exit_code and returns the same value
  typeset -i _ec
  printf "PrEx $(opa_date)    > %s\n" "$*"
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
  printf "PrEx $(opa_date) %3d< %s [%s]\n" "$_ec" "$*" "$_status"
  return $_ec
}

function opa_die {
  # log an error message and die
  typeset -i _ec="$1"
  shift 1
  opa_log 0 "$@"
  opa_exit $_ec
}

function opa_log {
  # Write a log messages ($2...) on standard output;
  # is $1 is lower than $opa_LOGMAX, the messages is written also onto $opa_LOGFILE
  if [[ $opa_LOGMAX -eq $opa_DISABLE_LOG ]] ; then
    return 0
  fi
  typeset -i _nlog="$1"
  shift 1
  for _message in "$@" ; do
    if ( $opa_BATCH ) && [[ $opa_LOGMAX -ne $opa_DISABLE_LOG_ON_LOGFILE && $_nlog -le $opa_LOGMAX ]] ; then
      printf '%s %02d (%s) : %s\n' "$(opa_date)" "$_nlog" "$opa_RUNID" "$_message" | tee -a $opa_LOGFILE
    else
      printf '%s %02d (%s) : %s\n' "$(opa_date)" "$_nlog" "$opa_RUNID" "$_message"
    fi
  done
}

function opa_exit {
  # Log message for script exit
  typeset -i _default_ec=${_opa_EXIT_CODE:-$?}
  typeset -i _ec=${1:-$_default_ec}
  #echo "opa_exit(${1}->${_default_ec}->${_ec})"
  typeset    _message=""
  if [[ "$1" == "" ]] ; then
    _message=" (from last executed command)"
  fi
  if $_opa_START_DONE && ! $_opa_EXIT_DONE ; then
    opa_log 0 "* EXIT(${opa_PROGNAME}) [$(opa_exitcode ${_ec})]$_message"
    _opa_EXIT_DONE=true
    _opa_EXIT_CODE=${_ec}
  fi
  exit $_ec
}


