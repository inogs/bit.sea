from __future__ import print_function
from sys import stderr
from subprocess import check_output

LOG_LEVELS = {1 : 'report', 2 : 'info', 3 : 'debug'}
COL_SPACE = 9

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
        self._lines = []
    
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
        return self.__verbose

    def separation_line(self):
        print('-'*self.__tty_cols)

    def set_verbosity_level(self, v):
        self.__verbose = int(v)

    def get_content(self):
        output = ""
        for l in self._lines:
            output += LOG_LEVELS[l[0]].upper() + ': ' + l[1] + '\n'
        return output

    def error(self, txt, *args, **kwargs):
        first_word = ('{:<' + str(COL_SPACE) +'}').format('ERROR:')
        txt = self._format_text(kwargs, first_word, txt) 
        print(txt, *args, file=stderr)


def add_log_level(l):
    def f(self, txt, *args, **kwargs):
        self._lines.append((l,txt.replace('\n', ' ')))
        first_word = ('{:<' + str(COL_SPACE) +'}').format(LOG_LEVELS[l].upper() + ':')
        if self.get_verbosity_level() >= l:
            txt = self._format_text(kwargs, first_word, txt)
            print(txt, *args, **kwargs)
    setattr(Log, LOG_LEVELS[i], f)

for i in LOG_LEVELS:
    add_log_level(i)
