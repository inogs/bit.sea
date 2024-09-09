Write your own harvester
---------------------------------
All the information reported in the previous section shall be enough to run the software using all its functionality. The following parts
are useful only if you want to develop your own harvester.


How the download script works
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The download script imports any python module that is inside the harvesters dir or that is imported by then. For example, if a file in
the harvesters directory imports another file in another directory, also that file will be imported by the downloader script. Then it
checks, in all the imported elements, if there are some class that ends with the world Harvester. These are the class that will be
executed. In other words, it does not matter the filesystem structure or the file names: you can have 20 harvester in a single file
or you can have some files in the harvesters directory that do not declare any harvester. The harvesters are simply the class that ends
their names with the string "Harvester". 

How the harvester should be implemented
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
You are completely free to implement the harvester in the way you prefer. All that you have to know is that the downloader script will
call a method called "harvest" on an object created from this class (or "rebuild" if in reset mode). This means that:

  - __init__ should not require any arguments (otherwise the downloader script would not be able to create an object for that harvester class)
  - two methods called "harvester" and "rebuild" shall be implemented in your harvester class. The first two arguments that will be passed
    to these methods are the general path where the file should be saved and a log object to send message to the user. Please read the following
    section about the log element.

Moreover, it is expected that every method (harvester and rebuild) returns a list of all the files that have been downloaded during the
execution of the method. While the return value of these methods is discard at the moment, in the future it could be used.

Do NOT use print
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Please, while implement your own harvester, do NOT use the print command. While in principle it would work, there are several disadvantages:

  - The message that you print on the screen with "print" will not be reported on the log file
  - For the very same reason, it will not be reported in the possible mail
  - The verbosity option will not affect the behaviour of print

For these reasons, it would be better to use the log object that the downloader script passes to the harvest and rebuild methods.
To use it, simplit call the method with the severity of the message (eg. log.info("My message") or log.debug("Another message")).

To know which level of severity are implemented in the log object, take a look at the Log class documentation.
