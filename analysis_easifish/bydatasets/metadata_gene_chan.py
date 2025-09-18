# all projection data genes and channels

def get_proj_metadata():
    """
    """
    proj_metadata = {
        'lt185': {'channels': ['r1v3_c0', 'r1v3_c2'], 
                  'proj_targets': ['RL', 'LM'], 
                 }, 
        'lt186': {'channels': ['r1v3_c0', 'r1v3_c2'], 
                  'proj_targets': ['RL', 'LM'],
                  'axis_polarity': [1,1,-1], # flip z or invert y axis
                 }, 
        'cdf03_c1-2_sub': {'channels': ['r1_c0', 'r1_c2'], 
                  'proj_targets': ['LM', 'RL'], 
                  'axis_polarity': [1,1,-1], # flip z or invert y axis
                 }, 
        'cdf04_c1-2_sub': {'channels': ['r1_c0', 'r1_c2'], 
                  'proj_targets': ['LM', 'RL'], 
                  'axis_polarity': [1,1,-1], # flip z or invert y axis
                 }, 
        
        
        
        'cdf03_c1-2': {'channels': ['r1_c0', 'r1_c1', 'r1_c2'], 
                  'proj_targets': ['LM', 'AM/PM', 'RL'], 
                  'axis_polarity': [1,1,-1], # flip z or invert y axis
                 }, 
        'cdf04_c1-2': {'channels': ['r1_c0', 'r1_c1', 'r1_c2'], 
                  'proj_targets': ['LM', 'AL', 'RL'], 
                  'axis_polarity': [1,1,-1], # flip z or invert y axis
                 }, 
        
        'cdfv0_c1-2': {'channels': ['r1_c0', 'r1_c1', 'r1_c4', ], 
                  'proj_targets': ['AL', 'AM/PM', 'LM', ], 
                 }, # RL in r2 but non-consistent z-step
        'cdfv1_c1-2': {'channels': ['r1_c0', 'r1_c1', 'r1_c4', 'r2_c1'], 
                  'proj_targets': ['AL', 'AM/PM', 'LM', 'RL'], 
                 }, # RL in r2; LM c4 bleeding issue; AL vCre (new - strong signals)    
    }
    
    
    default_target_colors = {
        'LM': 'red',
        'RL': 'blue',
        'AM/PM': 'green',
        'AL': 'orange',
    }
    
    for key in proj_metadata.keys():
        proj_metadata[key]['colors'] = [default_target_colors[target]
                                        for target in proj_metadata[key]['proj_targets']]
    
    return proj_metadata
