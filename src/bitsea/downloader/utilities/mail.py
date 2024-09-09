# Copyright (c) 2015 eXact Lab srl
# Author: Stefano Piani <stefano.piani@exact-lab.it>
from subprocess import Popen, PIPE
from time import time
from sys import version_info
from hashlib import md5
from base64 import b64encode

# Get the current python version
py_version = int(version_info[0])


# Unfortunately, python2 and python3 handle the strings
# in a different way. In python 2, a sequence of bytes
# and a string are the same thing, in python3 we need to
# perform a conversion
if py_version < 3:
    def str_to_bytes(txt):
        return txt
    def bytes_to_str(b):
        return b
else:
    def str_to_bytes(txt):
        return txt.encode('ASCII')
    def bytes_to_str(b):
        return b.decode()

# This is useful to send an attachment using email
def str_to_base64(txt):
    """Take a string and rewrite it in base64"""
    return bytes_to_str(b64encode(str_to_bytes(txt)))



def write_mail(sender, recipient_list, subject, text,
               attachment_name=None,
               attachment_content=None):

    """
    Send a mail to someone. If needed, a text file can be attached
    to the mail.
    This function relies on the sendmail utilities which should be
    properly configured.
    
    Args:
        - *sender*: The name of the sender as a string without spaces
        - *recipient_list*: A list of strings that are the addresses of 
          the recipients
        - *subject*: The subject of the mail as a string
        - *text*: The text of the mail as a string
        - *attachment_name* (optional): If it is not None, this is the
          name of the attached file that will be send with the mail.
          If this parameter is not none, also the attachment_content
          shall not be empty.
        - *attachment_content* (optional): The content of the attached
          file as a string. This content will be converted in Base64
          and will be attached to the mail.
    """
 
    recipient_str = ', '.join(recipient_list)
    m = md5()
    m.update(str_to_bytes(str(time())))
    boundary = m.hexdigest()
    
    # header of the message
    message_text  = 'From: {}\n'.format(sender)
    message_text += 'To: {}\n'.format(recipient_str)
    message_text += 'Subject: {}\n'.format(subject)
    message_text += 'Mime-Version: 1.0\n'
    message_text += 'Content-Type: multipart/mixed; boundary={}\n\n'.format(boundary)
    
    # the text
    message_text += '--{}\n'.format(boundary)
    message_text += 'Content-Type: text/plain; charset="US-ASCII"\n'
    message_text += 'Content-Transfer-Encoding: 7bit\n'
    message_text += 'Content-Disposition: inline\n\n'
    message_text += text
    
    # the attachment as a base64 string
    if attachment_name is not None:
        message_text += '\n\n'
        message_text += '--{}\n'.format(boundary)
        message_text += 'Content-Type: text/plain; '
        message_text += 'charset=US-ASCII; '
        message_text += 'name="{}"\n'.format(attachment_name)
        message_text += 'Content-Disposition: attachment; '
        message_text += 'filename="{}"\n'.format(attachment_name)
        message_text += 'Content-Transfer-Encoding: base64\n\n'

        message_text += str_to_base64(attachment_content)
        
    # End message
    message_text += '\n\n--{}--'.format(boundary)
    
    sendmail = Popen(['sendmail', '-t', '-oi'], stdin=PIPE, stdout=PIPE)
    sendmail.communicate(str_to_bytes(message_text))
