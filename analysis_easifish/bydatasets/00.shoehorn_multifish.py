"""
"""
"""Use python3
"""
import os
import shutil
import json
import glob


# ddir = "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/lt186_r2_autos1_flatfused.n5"
# ddir = "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/lt186_r4_autos1_flatfused.n5"
# ddir = "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/lt186_r3_autos1_flatfused.n5"
# ddir = "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/lt186_r6_autos1_flatfused.n5"
# ddir = "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/lt186_r7_autos1_flatfused.n5"
# ddir = "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/lt186_r1_autos1_flatfused.n5"

# ddir = "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/cdf03_c1-2_bino_r1_autos1_flatfused.n5"
# ddir = "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/cdf04_c1-2_bino_r1_autos1_flatfused.n5"

# ddir = "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/cdf03_c1-2_bino_r2_autos1_flatfused.n5"
# ddir = "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/cdf04_c1-2_bino_r2_autos1_flatfused.n5"

# ddir = "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/lt185_r2_autos1_flatfused.n5"
# ddir = "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/lt185_r3_autos1_flatfused.n5"
# ddir = "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/lt185_r4_autos1_flatfused.n5"
# ddir = "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/lt185_r5_autos1_flatfused.n5"
# ddir = "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/lt185_r7_autos1_flatfused.n5"
# ddir = "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/lt185_r6_autos1_flatfused.n5"
# ddir = "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/lt185_r1_autos1_flatfused.n5"

# ddir = "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/cdf03_c1-2_bino_r3_autos1_flatfused.n5"
# ddir = "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/cdf04_c1-2_bino_r3_autos1_flatfused.n5"
# ddir = "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/cdf03_c1-2_bino_r3_1_autos1_flatfused.n5"

# ddir = "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/lt172_r1_autos1_flatfused.n5"
# ddir = "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/lt172_r2_autos1_flatfused.n5"
# ddir = "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/lt172_r3_autos1_flatfused.n5"
# ddir = "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/lt172_r4_autos1_flatfused.n5"
# ddir = "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/lt172_r5_autos1_flatfused.n5"

# ddir = "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/test2_r1_autos1_flatfused.n5"

# ddirs = [
#     f"/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/273LU_aldh1l1Cre_Pcdhg_homo_autos1_flatfused_tile{i}.n5"
#     for i in range(12)
# ] + [
#     f"/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/296LU_c3OE_cortex_autos1_flatfused_tile{i}.n5"
#     for i in range(11)
# ] + [
#     f"/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/296LURU_pcdhg_homo_autos1_flatfused_tile{i}.n5"
#     for i in range(9)
# ] + [
#     f"/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/456RU_A1OE_cortex_autos1_flatfused_tile{i}.n5"
#     for i in range(5)
# ]

ddirs = [
    # "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/cdfv0_c1-2_r1_autos1_flatfused.n5",
    # "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/cdfv1_c1-2_r1_autos1_flatfused.n5",
    # "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/cdfv1_c1-2_r2_autos1_flatfused.n5",
    # "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/cdfv0_c2-2_r1_autos1_flatfused.n5",
    # "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/cdfv1_c2-1_r1_autos1_flatfused.n5",
    # "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/cdfv0_c1-1_r1_autos1_flatfused.n5",
    # "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/cdfv1_c1-1_r1_autos1_flatfused.n5",
    # "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/cdfv1_c1-3_r1_autos1_flatfused.n5",
    # "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/cdfv1_c2-2_r1_autos1_flatfused.n5",
    "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/cdfv1_c1-2_r1_REDO_autos1_flatfused.n5",
    "/u/home/f/f7xiesnm/project-zipursky/data/03_flatfused/cdfv0_c1-2_r1_REDO_autos1_flatfused.n5",
]

channels=['c0', 'c1', 'c2',] #'c3', 'c4']

for ddir in ddirs:
    print(f"->{ddir}")
    for chan in channels:
        print(f'-->{chan}')
        path = os.path.join(ddir, chan) 

        resolutions = [res for res in os.listdir(path) 
                        if res.startswith('s') and 
                        not os.path.islink(os.path.join(path, res))]

        resolution_integers = [int(res[1:]) for res in resolutions] 
        #### debug only !  resolution_integers = [int(res[1:-len('_bck')]) for res in resolutions] 

        new_resolutions = [f"s{ri+1}" for ri in resolution_integers]
        print(resolutions)
        print(resolution_integers)
        print(new_resolutions)

        ## rename dir, copy json, and modify json
        for ri, res in zip(resolution_integers, resolutions):
            src     = os.path.join(ddir, chan, res)
            dst     = os.path.join(ddir, chan, res+'_bck')
            jsn     = os.path.join(ddir, chan, res+'_bck', 'attributes.json')
            jsn_bck = os.path.join(ddir, chan, res+'_bck', 'attributes.json.bck')

            #### debug only ! jsn     = os.path.join(ddir, chan, res, 'attributes.json')
            #### debug only ! jsn_bck = os.path.join(ddir, chan, res, 'attributes.json.bck')

            ## rename dir
            if not os.path.isdir(dst):
                os.rename(src, dst)
            ## copy attributes
            if not os.path.isfile(jsn_bck):
                shutil.copy(jsn, jsn_bck)
            ## modify json
            with open(jsn_bck, 'r') as fh:
                a = json.load(fh)
                a['downsamplingFactors'] = [2**(ri+1),2**(ri+1),2**(ri)] 
                a['pixelResolution'] = [0.23, 0.23, 0.42] 
            with open(jsn, 'w') as fh:
                json.dump(a, fh)

        ## link to new resolution
        for res, newres in zip(resolutions, new_resolutions):
            dst = os.path.join(ddir, chan, res+'_bck')
            lnk = os.path.join(ddir, chan, newres)
            ## symlink
            if not os.path.isdir(lnk) and not os.path.islink(lnk):
                os.symlink(dst, lnk)

