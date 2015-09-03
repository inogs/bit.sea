The Downloader script
-------------------------

The downloader script is a software whose goal is to synchronize one or more remote folders on a common
path of the machine where the script is launched.

For every remote folder, there is an "harvester" which is in charge of download all the files of that particular
folder. The downloader script calls every harvester and then reports if any error occurred.

The downloader script expects that a common folder where all the data could be saved exists on the local machine.
Inside this common folder, every harvester can have its own subdirectory. There are two ways to communicate to the
downloader script the position of the common directory:

  - Using the -d flag in the command line interface
  - Changing the DEF_PATH variable at the beginning of the script and not setting the -d flag

The path of this directory will be passed to every harvester which is now in charge of save the remote files in a 
folder that take into account this path. While, in principle, an harvester is completely free of even ignore the
specified path, the usual behavior is to save the files in a subdirectory of the one passed with the -d flag; the
name of the subdirectory (or the subpath) is saved in a variable at the beginning of the file that define the
harvester. For example, for the BioFloatHarvester, the default subdir is "FLOAT_BIO".

Every harvester can run in two different ways: it can "harvest" or it can "rebuild". If the harvester is in the harvest
mode then it will simply update the archive: it will believe that on the local directory there is a copy of the remote dir
as it was some times before and it will try to update it at the current version. Please, keep in mind that if you edit or
remove some file inside the directory there is no guarantee that, calling the downloader script in the  harvester mode,
you will have a local copy of the remote folder. For example, let say that you have a file called "f1" written on April
and one called "f2" written on May. If you delete the file f1 it is possible that the harvester, looking at the f2 file,
believed that there is no need to download any file edited before May and therefore, it will not download again the f1 file.

Instead, if an harvester is run in "rebuild" mode, then it will download again every file from the remote dir, completely
discarding every edit on the local copy. At the end you will have a perfect synchronized copy of the remote dir, no matter
if the local dir was edited.


The command line interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Using the command line interface, it is possible to decide the behavior of the downloader script. Please, launch it with
the "-h" flag to see every possible option. Here you will find a more detailed description of what they do.

  - *-a* or *--alert*: It sends a mail to the address specified if there are some problems executing one of the
    harvester. Please have a look to the mail section to have a better understand of how this works
  - *-d* or *--directory*: It sets where is the local directory that should be synchronized
  - *-e* or *--execute-only*: It execute only the harvesters that contain the argument inside their name. For
    example, with the -e bio you will execute the harvester BioFloatHarvester but not the MyOceanHarvester. The
    argument is not case sensitive
  - *-g* or *--gid*: The program will be executed using the privileges associated with the group whose id is
    selected. If you want to know which groups your user belongs please use the command "groups <username>". If
    you then want to know the id of a particular group, look inside the file /etc/groups. All the files and the
    directories created by the script will belong to the group selected with this option
  - *-l* or *--log-file*: Set the name of the file that will be used to store a log of the execution. If empty,
    no file will be created and the log will not be saved (empty means -l "").
  - *-r* or *--report*: Exactly like the --alert flag, but the mail will be send even if the execution does not
    encounter any error. As usual, have a look at the mail section
  - *-u* or *--uid*: The program will be executed using the privileges associated with the user whose id is
    selected. If you want to know the id of a particular user, use the command "id -u <username>". Take into
    account that in a usual setup only root can change the user id of a process, so, if you are not logged as
    root, it is likely that this command will fail. If not, all the files and the directories created by the 
    script will belong to the user selected with this option
  - *-v* or *--verbose*: Set the verbosity of the script. These are the possible levels:
    
      + 0: Write something only if there are some errors
      + 1: Write the errors and when the execution of an harvester starts and ends
      + 2: Write the most common informations
      + 3: Write every single operation which is performed
    
    The default is 2. If you want more informations, please read the log section.
  - *--RESET*: To ensure that this flag is not activate by mistake, there is not short version of it. If set,
    the harvester will run in the rebuild mode instead that in the harvester mode.


The mail utility
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
With the command line, it is possible to specify an address (or a list of address split by a comma) of people
that should be warned if an execution encounter an error or, for testing purpose, even in every execution using
the flag -r.
THe script relies on the utility "sendmail". This utility is called by the command "sendmail" on the command line
interface. Please, ensure that this software is installed and properly configured before use the mail option.



