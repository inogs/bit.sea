#Helpers function to navigate the DOM

def get_subelements(node, tag):
    return [e for e in node.childNodes if (e.nodeType == e.ELEMENT_NODE and e.localName == tag)]

def get_node_attr(node, attribute):
    return node.attributes[attribute].value

