'''
Created on Mar 23, 2015

@author: ayan
'''
import functools
import warnings



class CannotFindPaddingError(Exception):
    
    base_message = 'The netCDF file appears to have conform to SGRID conventions, but padding values cannot be found.'
        
    def __str__(self):
        return self.base_message


class SGridNonCompliantError(Exception):
    
    base_message = 'This netCDF object does not appear to be SGRID compliant: {0}.'
    
    def __init__(self, dataset_obj):
        self.dataset_obj = dataset_obj
        
    def __str__(self):
        filepath = self.dataset_obj.filepath()
        error_message = self.base_message.format(filepath)
        return error_message
    
    
def deprecated(deprecated_function):
    @functools.wraps(deprecated_function)
    def new_func(*args, **kwargs):
        warnings.warn_explicit('Call to deprecated function: {0}'.format(deprecated_function.__name__),
                               category=DeprecationWarning,
                               filename=deprecated_function.__code__.co_filename,
                               lineno=deprecated_function.__code__.co_firstlineno + 1
                               )
        return deprecated_function(*args, **kwargs)
    return new_func
        