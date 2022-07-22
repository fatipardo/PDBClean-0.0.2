#
import os
#
def remove_file_defined_chain_from_list(list):
    """
    """
    remove_chid = []
    print("    Enter the file name containing the list of chain IDs you want removed from Standard Sequences.")
    user_input = input('File: ')
    if (os.path.isfile(user_input) == True):
        my_file = open(user_input)
        for line in my_file:
            remove_chid.append(line.strip())
    else:
        print("File does not exist.")
    list = remove_chid_from_list(list, remove_chid)
    return list

def remove_user_defined_chain_from_list(list):
    """
    """
    remove_chid = []
    user_input  = ""
    print("    Enter chain IDs of the chains you want removed. When done, enter DONE.")
    while (user_input != "DONE"):
        user_input = input('Chain ID: ')
        remove_chid.append(user_input)
    list = remove_chid_from_list(list, remove_chid)
    return list

def remove_chid_from_list(list, remove_list):
    """
    """
    for elt in remove_list:
        if elt in list:
            list.remove(elt)
#            del chlist[chid]
    return list

def show_list(list):
    """
    show_list
    """
    for elt in list:
        print(elt)

