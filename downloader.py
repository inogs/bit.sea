from __future__ import print_function

from sys import stderr, exit
from pkgutil import iter_modules
from importlib import import_module
from inspect import isclass, getmembers
from os import makedirs
from os.path import dirname, isdir, exists, realpath
from argparse import ArgumentParser, RawTextHelpFormatter

import traceback

import harvesters
from utilities.log_class import Log

# This is the default location where the database will be stored if
# that information is not specified when running the script
DEF_PATH = dirname(realpath(__file__)) +  '/data'



if __name__ == '__main__':

    # Read the command line arguments
    parser = ArgumentParser(formatter_class=RawTextHelpFormatter,
                            add_help=False)
    parser.add_argument('-a', '--alert', type=str, default="",
                        help='\nSend an email at the end of the '
                             'execution of the script in case of '
                             'error in one harvester\n\n')
    parser.add_argument('-d', '--directory', type=str, default=DEF_PATH,
                        help="Set the directory where the database " +
                             "will be stored\nDefault: "+DEF_PATH +"\n\n")
    parser.add_argument('-e', '--execute-only', type=str, default='',
                        help="Execute an harvester only if contains the " +
                             "argument in its name. This is particularly "
                             "useful for debug purposes.\nDefault:'' \n\n")
    parser.add_argument('-h', '--help', action='help',
                        help='\nShow this help and exit\n\n')
    parser.add_argument('-r', '--report', type=str, default="",
                        help='\nSend an email at the end of the '
                             'execution of the script with the '
                             'log of all the operations\n\n')
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

    # Create a log object
    log = Log(verbose_level)

    # Check if the directory of the database exists. If not, create it
    if exists(database_path):
        if isdir(database_path):
           log.info(database_path + " exists and it is a directory!")
        else:
            log.error(database_path + " exists but it is not a"
                      "directory!")
            exit(1)
    else:
        if reset_mode == False:
            log.queerness("Reset mode is not active, but there is no "
                          "archive to update! The script will work "
                          "in reset mode even if it is not specified "
                          "in the command line interface")
            reset_mode = True
        makedirs(database_path)
        log.info("Directory " + database_path + " created")

    # Look every module inside the directory of the harvesters
    list_of_harvester_modules = []
    harvesters_path = harvesters.__path__
    all_modules = iter_modules(harvesters_path)
    for _, modname, ispkg in all_modules:
        if not ispkg:
            mod_path = 'harvesters.' + modname
            current_mod = import_module(mod_path)
            log.info('Imported module ' + mod_path)
            list_of_harvester_modules.append(current_mod)

    # Look for every class with a name that ends with "Harvester" in that
    # modules
    list_of_harvester_classes = []
    for current_module in list_of_harvester_modules:
        for obj_name, obj in getmembers(current_module):
            if isclass(obj) and obj_name[-9:]=="Harvester":
                list_of_harvester_classes.append((obj_name, obj))


    # Now, for every class, we call the method harvest (or rebuild if we
    # are running in reset mode. In case of an error, it will be reported
    # at the end of the process

    errors = []
    for harvester_name, Harvester in list_of_harvester_classes:
        if execution_filter in harvester_name:
            log.separation_line()
            log.achievement('Running ' + harvester_name)
            harvester = Harvester()
            try:
                if reset_mode:
                    harvester.rebuild(database_path, log)
                else:
                    harvester.harvest(database_path, log)
                log.achievement(harvester_name + ' sucessfully executed!')            
            except:
                exception = traceback.format_exc()
                errors.append((harvester_name, exception))
                log.error("An error occurred executing the "
                          "harvester " + harvester_name + "!")
    
    if len(errors) != 0:
        log.separation_line()
        for harvester, exception in errors:
            log.error("Error executing " + harvester +
                      ":\n" + exception, split_lines=False)
        exit(1)
