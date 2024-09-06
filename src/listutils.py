#
import os
#
def remove_file_defined_chain_from_list(list):
    """
    Removes specified chain IDs from list based on user's file input.

    The user is prompted to enter the file name containing the chain IDs they want to remove.
    The chain IDs in the file will be removed from the chain ID list.

    Parameters:
    -----------
    list : list
        Contains all the chain IDs from CIF(s)

    Returns:
    --------
    list : list
        Updated list without the chain ID from user's input
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
    Removes the chain ID from list based off of user's input of chain ID

    The user is prompted to input the chain ID which they wish to remove

    Parameter:
    ----------
    list : list
        Contains all the chain IDs from CIF(s)

    Returns:
    --------
    list : list
        Updated list without the chain ID from user's input
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
    Removes the chain ID from the list

    Parameters:
    -----------
    list : list
        contains the chain IDs
    remove_list : list
        contains the chain IDs to be removed

    Returns:
    --------
    list : list
        Updated list without the chain ID from user's input
    """
    for elt in remove_list:
        if elt in list:
            list.remove(elt)
    return list

def show_list(list):
    """
    This function prints each item in a list.

    Parameters:
    -----------
    list : list
    	A list of items to be printed.

    Returns:
    --------
    None
    """
    for elt in list:
        print(elt)

