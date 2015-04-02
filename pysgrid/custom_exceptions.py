'''
Created on Mar 23, 2015

@author: ayan
'''


class CannotFindPadding(Exception):
    
    base_message = 'The netCDF file appears to have conform to SGRID conventions, but padding values cannot be found.'
        
    def __str__(self):
        error_message = self.base_message.format(self.padding_str)
        return error_message


class SGridNonCompliant(Exception):
    
    base_message = 'This netCDF object does not appear to be SGRID compliant: {0}.'
    
    def __init__(self, dataset_obj):
        self.dataset_obj = dataset_obj
        
    def __str__(self):
        error_message = self.base_message.format(self.dataset_obj)
        return error_message
    
    
class DimensionMismatch(Exception):
    
    def __init__(self, *args):
        self.args = args
        
    def __str__(self):
        error_message = ('The is a dimension mismatch between arrays with shapes {0}. '
                         'Arrays must have the same same to use this function.')
        array_shapes = [arr.shape for arr in self.args]
        array_shapes_count = len(array_shapes)
        if array_shapes_count == 2:
            shape_str = 'and '.join(array_shapes)
        else:
            shape_str_intermediate = ', '.join(array_shapes).rsplit(',', 1)
            shape_str = ', and'.join(shape_str_intermediate)
        filled_message = error_message.format(shape_str)
        return filled_message
        