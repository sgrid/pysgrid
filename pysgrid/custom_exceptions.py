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
    
    base_message = 'This netCDF object derived from the dataset at {0} does not appear to be SGRID compliant.'
    
    def __init__(self, dataset_path):
        self.dataset_path = dataset_path

    def __str__(self):
        error_message = self.base_message.format(self.dataset_path)
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
        