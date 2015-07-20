from __future__ import print_function
from sys import stderr
from subprocess import check_output

class Log(object):
    def __init__(self, verbose_level=1):
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
    
    def __restruct_in_lines(self, first_word, txt):
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

    def get_verbosity_level(self):
        return self.__verbose

    def separation_line(self):
        print('-'*self.__tty_cols)

    def set_verbosity_level(self, v):
        self.__verbose = int(v)

    def info(self, txt, *args, **kwargs):
        txt = self.__restruct_in_lines('INFO:    ', txt)
        if self.__verbose > 2:
            print(txt, *args, **kwargs)

    def queerness(self, txt, *args, **kwargs):
        if 'split_lines' in kwargs:
            if kwargs['split_lines'] == True:
                txt = self.__restruct_in_lines('STRANGE: ', txt)
            else:
                txt = '\n'.join([' '*9 + l for l in txt.split('\n')])
                txt = 'STRANGE:' + txt[8:]
            del kwargs['split_lines']
        else:
            txt = self.__restruct_in_lines('STRANGE: ', txt)

        if self.__verbose > 1:
            print(txt, *args, **kwargs)

    def achievement(self, txt, *args, **kwargs):
        txt = self.__restruct_in_lines('REPORT:  ', txt)
        if self.__verbose > 0:
            print(txt, *args, **kwargs)

    def error(self, txt, *args, **kwargs):
        if 'split_lines' in kwargs:
            if kwargs['split_lines'] == True:
                txt = self.__restruct_in_lines('ERROR:   ', txt)
            else:
                txt = '\n'.join([' '*9 + l for l in txt.split('\n')])
                txt = 'ERROR:' + txt[6:]
            del kwargs['split_lines']
        else:
            txt = self.__restruct_in_lines('ERROR:   ', txt)
            
        print(txt, *args, file=stderr)
