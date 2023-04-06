""" Functions that are used by multiple scripts bundled """

def flatten(l_l):
    """
    flattens a list of lists consisting of strings 

    :param l_l: list of lists
    :type values: list
    :return: flattened list
    :rtype: list
    """
    flattened_list = []
    for i in l_l:
        if isinstance(i,list): 
            flattened_list.extend(flatten(i))
        else: 
            flattened_list.append(i)
    return flattened_list