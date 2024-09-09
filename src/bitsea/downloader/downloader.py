# Copyright (c) 2015 eXact Lab srl
# Author: Stefano Piani <stefano.piani@exact-lab.it>
from __future__ import print_function

from sys import stderr, exit
from pkgutil import iter_modules
from importlib import import_module
from inspect import isclass, getmembers
from os import makedirs, setuid, setgid
from os.path import dirname, isdir, exists, realpath
from argparse import ArgumentParser, RawTextHelpFormatter

import traceback

import harvesters
from utilities.log_class import Log
from utilities.mail import write_mail

# This is the default location where the database will be stored if
# that information is not specified when running the script
#DEF_PATH = dirname(realpath(__file__)) +  '/data'
DEF_PATH = "/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE"
LOG_FILE = "log/downloader.log"

if __name__ == '__main__':

    # Read the command line arguments
    parser = ArgumentParser(formatter_class=RawTextHelpFormatter,
                            add_help=False)
    parser.add_argument('-a', '--alert', type=str, default="",
                        help='Send an email at the end of the '
                             'execution of the script in case of '
                             'error in one harvester\n\n')
    parser.add_argument('-d', '--directory', type=str, default=DEF_PATH,
                        help="Set the directory where the database " +
                             "will be stored\nDefault: "+DEF_PATH +"\n\n")
    parser.add_argument('-e', '--execute-only', type=str, default='',
                        help="Execute an harvester only if contains the " +
                             "argument in its name. This is particularly "
                             "useful for debug purposes.\nDefault:'' \n\n")
    parser.add_argument('-g', '--gid', type=int, default=-1,
                        help="\nSet the group that will execute this script " +
                             "by its id.\nPlease, be sure to have the " +
                             "permission to do this kind of operation\n\n")
    parser.add_argument('-l', '--log-file', type=str, default=LOG_FILE,
                        help="The relative or absolute path of the file " +
                             "where the log will be saved. If the string " +
                             "is empty, no file will be written\n" +
                             "Default: " + str(LOG_FILE) + "\n\n")
    parser.add_argument('-h', '--help', action='help',
                        help='\nShow this help and exit\n\n')
    parser.add_argument('-r', '--report', type=str, default="",
                        help='Send an email at the end of the '
                             'execution of the script with the '
                             'log of all the operations\n\n')
    parser.add_argument('-u', '--uid', type=int, default=-1,
                        help="\nRun this script as the user with the uid set " +
                             "in this option.\nPlease, be sure to have the " +
                             "permission to do this kind of operation\n\n")
    parser.add_argument('-v', '--verbosity', type=int, default=2,
                        help="Define the level of verbosity\n"
                             "Default: 2\n\n")
    parser.add_argument('--RESET', action="store_true",
                        help="\nDelete all the archive and rebuild it "+
                             "from scratch\n\n")

    argv = parser.parse_args()
    database_path = argv.directory
    verbose_level = argv.verbosity
    reset_mode = argv.RESET
    execution_filter = argv.execute_only

    if argv.report != "":
        report_to = argv.report.split(',')
    else:
        report_to = []

    if argv.alert != "":
        alert_to = argv.alert.split(',')
    else:
        alert_to = []

    if argv.log_file == "":
        log_file = None
    else:
        log_file = argv.log_file
    
    # Create a log object
    log = Log(verbose_level, log_file)

    # Try to set gui and uid
    if argv.gid != -1:
        try:
            setgid(argv.gid)
        except:
            exception = traceback.format_exc()
            log.info('Impossible change the group id of the program. '
                     'The execution will continue without changing it.')
            log.debug("This is the traceback:\n" + exception,
                                            split_lines = False)
    if argv.uid != -1:
        try:
            setuid(argv.uid)
        except:
            exception = traceback.format_exc()
            log.info('Impossible change the user id of the program. '
                     'The execution will continue without changing it.')
            log.debug("This is the traceback:\n" + exception,
                                            split_lines = False)
            

    # Check if the directory of the database exists. If not, create it
    if exists(database_path):
        if isdir(database_path):
           log.debug(database_path + " exists and it is a directory!")
        else:
            log.error(database_path + " exists but it is not a"
                      "directory!")
            exit(1)
    else:
        if reset_mode == False:
            log.info("Reset mode is not active, but there is no "
                          "archive to update! The script will work "
                          "in reset mode even if it is not specified "
                          "in the command line interface")
            reset_mode = True
        makedirs(database_path)
        log.debug("Directory " + database_path + " created")

    # Look every module inside the directory of the harvesters
    list_of_harvester_modules = []
    harvesters_path = harvesters.__path__
    all_modules = iter_modules(harvesters_path)
    for _, modname, ispkg in all_modules:
        if not ispkg:
            mod_path = 'harvesters.' + modname
            current_mod = import_module(mod_path)
            log.debug('Imported module ' + mod_path)
            list_of_harvester_modules.append(current_mod)

    # Look for every class with a name that ends with "Harvester" in that
    # modules
    list_of_harvester_classes = []
    for current_module in list_of_harvester_modules:
        for obj_name, obj in getmembers(current_module):
            if isclass(obj) and obj_name[-9:]=="Harvester":
                list_of_harvester_classes.append((obj_name, obj))


    # Take only the harvesters that contains the filter parameter in ther
    # names 
    harvester_classes_filtered = [h for h in list_of_harvester_classes if 
                                  execution_filter.lower() in h[0].lower()]
    if len(harvester_classes_filtered) == 0:
        log.info('No harvester to be executed! Program will exit!')

    # Now, for every class, we call the method harvest (or rebuild if we
    # are running in reset mode. In case of an error, it will be reported
    # at the end of the process
    errors = []
    for harvester_name, Harvester in harvester_classes_filtered:
        log.separation_line()
        log.report('Running ' + harvester_name)
        harvester = Harvester()
        try:
            if reset_mode:
                harvester.rebuild(database_path, log)
            else:
                harvester.harvest(database_path, log)
            log.report(harvester_name + ' sucessfully executed!')            
        except:
            exception = traceback.format_exc()
            errors.append((harvester_name, exception))
            log.error("An error occurred executing the "
                      "harvester " + harvester_name + "!")
            log.debug("This is the traceback:\n" + exception,
                                            split_lines = False)

    # Report the errors
    if len(errors) != 0:
        log.separation_line()
        for harvester, exception in errors:
            log.error("Error executing " + harvester +
                      ":\n" + exception, split_lines = False)
        recipients = report_to + alert_to
        if len(recipients) > 0:
            message_text  = "Something went wrong executing the "
            message_text += "downloader script. Please, see the "
            message_text += "attached log!"
            write_mail("Downloader", recipients, "Something went wrong!!!",
                   message_text, LOG_FILE, log.get_content())
        exit(1)
    else:
        # Send the mail to the people in the report list
        if len(report_to) > 0:
            message_text  = "Please, find attached a complete report "
            message_text += "of the last execution of the downloader "
            message_text += "script"
            write_mail("Downloader", report_to, "Downloader report",
                   message_text, LOG_FILE, log.get_content())
