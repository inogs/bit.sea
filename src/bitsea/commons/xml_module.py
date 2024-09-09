# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
#Helpers function to navigate the DOM

def get_subelements(node, tag):
    return [e for e in node.childNodes if (e.nodeType == e.ELEMENT_NODE and e.localName == tag)]

def get_node_attr(node, attribute):
    output = None
    try:
        output = node.attributes[attribute].value
    except KeyError:
        output = None
    return output

