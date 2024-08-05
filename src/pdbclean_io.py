import sys, os, shutil, datetime
#

def check_project(projdir = None, level = 'top', action = 'create', verbose = True):
    """
    check_project
    """
    if projdir is None:
        print("Please provide a project directory path")
    else:
        dirname = projdir
        if(level!='top'):
            dirname = dirname+'/'+level
        if(action == 'create'):
            create_dir(dirname, verbose=verbose)
        elif(action == 'clean'):
            clean_dir(dirname, verbose=verbose)
        elif(action == 'delete'):
            delete_dir(dirname, verbose=verbose)

def create_dir(dirpath, verbose=True):
    """
    """
    if not os.path.exists(dirpath):
        os.mkdir(dirpath)
        if verbose:
            now = datetime.datetime.now()
            f = open(dirpath+'/info.txt', 'w')
            f.write('directory created on {0}'.format(now))
            f.close()
    else:
        if verbose:
            print('{0} already exists, with content:'.format(dirpath))
            print(os.listdir(dirpath))
            
def clean_dir(dirpath, verbose=True):
    """
    """
    if os.path.exists(dirpath):
        listfile = (file for file in os.listdir(dirpath) if os.path.isfile(os.path.join(dirpath, file)))
        if verbose:
            print('Cleaning {0}...'.format(dirpath))
        for f in listfile:
            os.remove(dirpath+'/'+f)

def delete_dir(dirpath, verbose=True):
    """
    """
    if os.path.exists(dirpath):
        shutil.rmtree(dirpath)
        if verbose:
            print('Deleting {0}...'.format(dirpath))
