merfish_datasets = {
    ## batch 1
    # exp 1-1
    'P14NRa_ant': "merfish_06142023/ant/region0",
    'P28NRa_ant': "merfish_06142023/ant/region1",
    
    # exp 1-2
    'P14NRa_pos': "merfish_06142023/pos/region0",
    'P28NRa_pos': "merfish_06142023/pos/region1",

    ## batch 2
    # exp 2-1
    'P21NRa_ant': "merfish_20231114/region0",
    'P21DRa_ant': "merfish_20231114/region2",
    'P28DRa_ant': "merfish_20231114/region1",
    
    # exp 2-2
    'P21NRa_pos': "merfish_20231120/region0",
    'P21DRa_pos': "merfish_20231120/region1",
    'P28DRa_pos': "merfish_20231120/region2",
    
    ## batch 3
    # exp 3-1
    'P14NRb_ant': "merfish_202404051211/region_0", 
    'P28NRb_ant': "merfish_202404051211/region_1",
    
    # exp 3-2
    'P14NRb_pos': "merfish_202404051214/region_0", 
    'P28NRb_pos': "merfish_202404051214/region_1",
    
    # exp 4-1
    'P21NRb_ant': "merfish_202404091526/region_2",
    'P21DRb_ant': "merfish_202404091526/region_1",
    'P28DRb_ant': "merfish_202404091526/region_0",
    
    # exp 4-2
    'P21NRb_pos': "merfish_202404091601/region_0",
    'P21DRb_pos': "merfish_202404091601/region_2",
    'P28DRb_pos': "merfish_202404091601/region_1",    
    
    ## batch 4
    # exp merfish_202503111431
    'P21NRc_ant': "merfish_202503111431/region_AP21",  # V1 good
    
    'P8NRa_ant':  "merfish_202503111431/region_AP81", 
    'P8NRb_ant':  "merfish_202503111431/region_AP82", 
    'P8NRc_ant':  "merfish_202503111431/region_AP83", 
    'P8NRd_ant':  "merfish_202503111431/region_AP84", 
    
    # exp merfish_202503111242
    'P21NRc_pos': "merfish_202503111242/region_PP21",   # V1 L2/3 cut off (imaging glitcsh)
    'R5_pos':     "merfish_202503111242/region_R5",     # V1 L2/3 cut off (imaging glitch)
    
    'P8NRa_pos':  "merfish_202503111242/region_PP81", 
    'P8NRb_pos':  "merfish_202503111242/region_PP82", 
    'P8NRc_pos':  "merfish_202503111242/region_PP83", 

    # exp merfish_202503191505 (repeated experiments) same animals diff. sections? 
    'P21NRc_ant2': "merfish_202503191505/region_P21R",    # V1 good
    
    'P8NRa_ant2':  "merfish_202503191505/region_P81", 
    'P8NRb_ant2':  "merfish_202503191505/region_P82", 
    'P8NRc_ant2':  "merfish_202503191505/region_P83", 
    'P8NRd_ant2':  "merfish_202503191505/region_P84", 
    
    # exp merfish_202503191540
    'P21NRc_pos2': "merfish_202503191540/region_P21p1",   # V1 good
    'P21NRc_pos2a': "merfish_202503191540/region_P21P1a", # V1 L2/3 cut off (sample glitch)

    'P8NRa_pos2':  "merfish_202503191540/region_P8p1", 
    'P8NRb_pos2':  "merfish_202503191540/region_P8p2", 
    'P8NRc_pos2':  "merfish_202503191540/region_P8p3", 
    'P8NRd_pos2':  "merfish_202503191540/region_P8p4", 
}

merfish_datasets['P14NR_ant'] = merfish_datasets['P14NRa_ant']
merfish_datasets['P28NR_ant'] = merfish_datasets['P28NRa_ant']

merfish_datasets['P14NR_pos'] = merfish_datasets['P14NRa_pos']
merfish_datasets['P28NR_pos'] = merfish_datasets['P28NRa_pos']

merfish_datasets['P21NR_ant'] = merfish_datasets['P21NRa_ant']
merfish_datasets['P21DR_ant'] = merfish_datasets['P21DRa_ant']
merfish_datasets['P28DR_ant'] = merfish_datasets['P28DRa_ant']

merfish_datasets['P21NR_pos'] = merfish_datasets['P21NRa_pos']
merfish_datasets['P21DR_pos'] = merfish_datasets['P21DRa_pos']
merfish_datasets['P28DR_pos'] = merfish_datasets['P28DRa_pos']


# note that the index and the labels are pretty complicated but that is true
merfish_datasets_params = { 
    'P14NRa_ant': {'rotation': 90, 'flip': False, 'rotation2': 45},
    'P14NRa_pos': {'rotation': 90, 'flip': False, 'rotation2': 45},

    'P14NRb_ant': {'rotation': 90, 'flip': False, 'rotation2': 45},
    'P14NRb_pos': {'rotation': 90, 'flip': False, 'rotation2': 45},
    
    
    # right hemi + 45; left hemi -45
    'P28NRa_ant': {'rotation':  90, 'rotation2': -45},
    'P28NRa_pos': {'rotation': -90, 'rotation2':  45},
    
    'P28DRa_ant': {'rotation': -90, 'rotation2': -45},
    'P28DRa_pos': {'rotation': -90, 'rotation2': -45},    
    
    'P28NRb_ant': {'rotation':  90, 'rotation2': -45},
    'P28NRb_pos': {'rotation':  90, 'rotation2': -45},
    
    'P28DRb_ant': {'rotation': -30, 'rotation2':  60},
    'P28DRb_pos': {'rotation':-200, 'rotation2':  60},
    
    
    
    # P8
    'P8NRa_ant2': {'rotation':  180, 'flip': False, 'rotation2': 45},
    'P8NRb_ant2': {'rotation': -160, 'flip':  True, 'rotation2': 45},
    'P8NRc_ant2': {'rotation':  180, 'flip':  True, 'rotation2': 45},
    'P8NRd_ant2': {'rotation':  180, 'flip':  True, 'rotation2': 45},
    
    'P8NRa_pos2': {'rotation':  180, 'flip': False, 'rotation2': 45},
    'P8NRb_pos2': {'rotation':   10, 'flip': False, 'rotation2': 45},
    'P8NRc_pos2': {'rotation':  180, 'flip': False, 'rotation2': 45},
    'P8NRd_pos2': {'rotation':  160, 'flip': False, 'rotation2': 45},
    
     # P21 needs update
    'P21NRc_ant': {'rotation':  180, 'flip': False, 'rotation2': 45}, 
    'P21NRc_pos': {'rotation':  180, 'flip': False, 'rotation2': 45}, 
    'R5_pos':     {'rotation':  180, 'flip': False, 'rotation2': 45}, 

    'P21NRc_ant2':  {'rotation':  180, 'flip': False, 'rotation2': 45}, 
    'P21NRc_pos2':  {'rotation':  180, 'flip': False, 'rotation2': 45},
    'P21NRc_pos2a': {'rotation':  180, 'flip': False, 'rotation2': 45},
    
    # previous P21
    'P21NRa_ant': {'rotation':  180, 'flip': False, 'rotation2': 45}, 
    'P21NRa_pos': {'rotation':  180, 'flip': False, 'rotation2': 45}, 
    'P21NRb_ant': {'rotation':  180, 'flip': False, 'rotation2': 45}, 
    'P21NRb_pos': {'rotation':  180, 'flip': False, 'rotation2': 45}, 
    
    # P21DR
    'P21DRa_ant': {'rotation':  180, 'flip': False, 'rotation2': 45}, 
    'P21DRa_pos': {'rotation':  180, 'flip': False, 'rotation2': 45}, 
    'P21DRb_ant': {'rotation':  180, 'flip': False, 'rotation2': 45}, 
    'P21DRb_pos': {'rotation':  180, 'flip': False, 'rotation2': 45}, 

}

ensemble_p21nr = [
    'P21NRb_ant',
    'P21NRb_pos',
    
    'P21NRc_ant',
    'P21NRc_ant2',
    'P21NRc_pos2',
]

ensemble_p21dr = [
    'P21DRa_ant',
    'P21DRa_pos',
    
    'P21DRb_ant',
    'P21DRb_pos',
]