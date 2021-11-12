
def get_database():
    import cv_analyzer, os
    datadir = cv_analyzer.__path__[0]+'/database/'
    files = os.listdir(datadir)
    files.sort()
    files = [os.path.join(datadir, i) for i in files if 'json' in i]
    return files
