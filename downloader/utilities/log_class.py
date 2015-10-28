# Copyright (c) 2015 eXact Lab srl
# Author: Stefano Piani <stefano.piani@exact-lab.it>
from __future__ import print_function
from sys import stderr
from subprocess import check_output
from datetime import datetime
from os.path import exists

# By definition, the -1 level is error and this can not
# be changed. The other levels are definited here
LOG_LEVELS = {1 : 'report', 2 : 'info', 3 : 'debug'}

# This is the space between the type of the message 
# (for example "info") and the content of the message
COL_SPACE = 9


class Log(object):
    """
    Log (short for Logger) is a class whose object are in charge of
    manage the output that shall be printed on the standard output and
    standard error, display it in a proper way, save it in memory or,
    optionally, also in a file.
    
    A log object store a verbose level that influence only what kind
    of information is printed on the standard output. In any case, the
    log will store in memory all the informations. The verbose level is
    an integer; a message will be print on the standard output only if
    its priority is greater or equal to the verbosity level of the log. 
    
    The possible levels of the log are stored in the LOG_LEVELS entries:
    it associates to an integer number a string which is the name of the
    level. That string is the name of the method of the log objects.
    
    Another level of severity it is also implemented even if it is not
    included in the LOG_LEVELS dictionary. It is the level "Error". This
    it is a little bit different from the others level for two reasons:
    
      - It will be printed no matter of the verbosity level
      - It will be printed on the standard error
    
    Args:
        - *verbose_level*: An integer that set the verbose level of the
          logger object
        - *file_log*: The name of the file where the log will be stored

    Returns:
        - *log*: a logger object 
    """ 
    def __init__(self, verbose_level, file_log=None):
        self.__verbose = verbose_level
        
        # This is a not portable way to know the dimension
        # of the terminal calling the stty utility
        try: 
            stty = check_output(['stty', 'size']).split()
            cols = int(stty[1])
            if cols > 100: cols=100
        except:
            cols = 100
        self.__tty_cols = cols

        # A memory buffer to store the log lines
        self._lines = []
        
        # A file where the log will be saved line by line
        self._file_log = file_log
        
        # Write the first line of the log
        if exists(file_log):
            skip_lines = True
        else:
            skip_lines = False
        if file_log is not None:
            with open(file_log, 'a') as f:
                if skip_lines:
                    f.write('\n\n')
                f.write('*' * COL_SPACE)
                f.write(' Starting a new execution at ')
                f.write(datetime.now().strftime("%I:%M%p on %B %d %Y"))
                f.write("\n")
    
    def _write_to_log(self, txt):
        # Save a string on the log file
        if self._file_log is not None:
            with open(self._file_log, 'a') as f:
                f.write(txt)
                f.write('\n')
    
    def _restruct_in_lines(self, first_word, txt):
        # Given a word and a some text, return a text with
        # the first world on the left and the text in another
        # coulumn on the right
        first_col = len(first_word)
        second_col = self.__tty_cols - first_col
        old_lines = txt.split('\n')
        new_lines = ['']
        for l in old_lines:
            words = l.split()
            for word in words:
                if len(new_lines[-1]) + len(word) + 1 < second_col:
                    new_lines[-1] = new_lines[-1] + ' ' + word
                else:
                    new_lines.append(' ' + word)
            new_lines.append('')
        # The following command removes the last empty line which
        # has been just added
        new_lines = new_lines[:-1]
        
        # Add empty space at the beginning. Moreover, we remove the
        # first char of every line which is, by construction, an
        # undesired whitespace
        for i in range(len(new_lines)):
            new_lines[i] = ' '*first_col + new_lines[i][1:]

        new_lines[0] = first_word + new_lines[0][first_col:]
        return '\n'.join(new_lines)

    def _format_text(self, kwargs, first_word, txt):
        if 'split_lines' in kwargs:
            if kwargs['split_lines'] == True:
                txt = self._restruct_in_lines(first_word, txt)
            else:
                txt = '\n'.join([' '*9 + l for l in txt.split('\n')])
                txt = first_word + txt[9:]
            del kwargs['split_lines']
        else:
                txt = self._restruct_in_lines(first_word, txt)            
        return txt

    def get_verbosity_level(self):
        '''Get the current level of verbosity of the logger'''
        return self.__verbose

    def separation_line(self):
        '''Print a line on screen to separate two sections'''
        if self.__verbose>=1:
            print('-'*self.__tty_cols)

    def set_verbosity_level(self, v):
        '''Change the verbosity level of the logger'''
        self.__verbose = int(v)

    def get_content(self, level=0):
        '''
        Return a string that contains all the information that
        the log has generated. Every message of the log is
        reported in one line. With the optional argument level
        is possible to select the verbosity of the content
        '''
        
        output = ""
        for l in self._lines:
            if l[0] == -1:
                output += 'ERROR: ' + l[1] + '\n'
            elif l[0] > level:
                output += LOG_LEVELS[l[0]].upper() + ': ' + l[1] + '\n'
        return output

    def error(self, txt, *args, **kwargs):
        '''
        Write the passed string on the standard error, putting an "ERROR:"
        string before. Moreover, if the string is too long, the text will
        be reformatted to be suitable for the console of the user. If you
        do not want this (for example because you are printing a traceback)
        you can disable this behaviour with the optional argument.
        split_line=False. All the other optional arguments will be passed
        directly to the print function that will be invoked by this method.
        
        Args:
            - *txt*: an error message
        
        '''

        # Save the message in the local buffer
        self._lines.append((-1, txt.replace('\n', ' ')))

        # Save the message in the log file
        self._write_to_log("ERROR: " + txt.replace('\n', ' '))
        
        # Print the message on the screen
        first_word = ('{:<' + str(COL_SPACE) +'}').format('ERROR:')
        txt = self._format_text(kwargs, first_word, txt) 
        print(txt, *args, file=stderr)



# This is a little bit complex. The function add_log_level
# add a new method for the class Log. The method added is the
# f function which is defined inside the very same function.
# l is the level of the log
def add_log_level(l):
    def f(self, txt, *args, **kwargs):
        # Save the message in the local buffer
        self._lines.append((l,txt.replace('\n', ' ')))

        # Save the message in the log file
        self._write_to_log(LOG_LEVELS[l].upper() + ": " + txt.replace('\n', ' '))

        first_word = ('{:<' + str(COL_SPACE) +'}').format(LOG_LEVELS[l].upper() + ':')
        if self.get_verbosity_level() >= l:
            txt = self._format_text(kwargs, first_word, txt)
            print(txt, *args, **kwargs)
    setattr(Log, LOG_LEVELS[l], f)


# Now we add a new method for every element of the LOG_LEVELS
# dictionary
for i in LOG_LEVELS:
    add_log_level(i)
