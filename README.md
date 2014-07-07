figs_base has a bunch of functions that can be used to create figures suitable for IEEE and ACM papers. 

Easiest to import as 

    import figs_base as fg
    
Or you might want to write a figs.py with your figure code. In this, use

   from figs_base import *

This is easier as several of the figs_base functions are intended to be called from other functions which compute specific figures for specific needs...

stats has some basic statistics functions.