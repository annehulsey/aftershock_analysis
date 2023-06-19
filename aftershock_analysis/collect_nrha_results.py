from .base import *


def store_building_geometry(results_filename, hf_group, n_stories, n_bays, story_heights, bay_widths):
    ## Creates arrays with the geometry meta data of the frame in the group with name hf_group inside the h5
    # file results_filename
    # INPUTS
    #       results_filename: name of h5 file to store frame information
    #       hf_group: name of the group to store builidng metadata
    #       n_stories: number of stories in frame
    #       n_bays: number of bays in frame
    #       story_heights: np.array with the story heights
    #       bay_widths: np.array with the bay lengths

    with h5py.File(results_filename, 'r+') as hf:
        hf_group = hf[hf_group]

        hf_group.attrs['n_stories'] = n_stories
        hf_group.attrs['n_floors'] = n_stories + 1

        hf_group.attrs['n_bays'] = n_bays
        hf_group.attrs['n_columns'] = n_bays + 1

        bay_widths = bay_widths * 12
        story_heights = story_heights * 12

        hf_group.attrs['bay_widths'] = bay_widths
        hf_group.attrs['story_heights'] = story_heights

        hf_group.attrs['building_width'] = np.sum(bay_widths)
        hf_group.attrs['building_height'] = np.sum(story_heights)

        # store the original geometry of each column
        columns = np.zeros(((n_bays + 1) * n_stories, 2, 2))
        i_element = 0
        for i_story in range(n_stories):
            for i_beam in range(n_bays + 1):
                # x values of columns
                columns[i_element, :, 0] = np.sum(bay_widths[:i_beam])
                for i_end in range(2):
                    # y values of columns
                    columns[i_element, i_end, 1] = np.sum(story_heights[:i_story + i_end])
                i_element = i_element + 1
        key = 'column_geometry'
        hf_group.create_dataset(key, data=columns)

        # store the original geometry of each beam
        beams = np.zeros((n_bays * n_stories, 2, 2))
        i_element = 0
        for i_story in range(n_stories):
            for i_beam in range(n_bays):
                # y values of beams
                beams[i_element, :, 1] = np.sum(story_heights[:i_story+1])
                for i_end in range(2):
                    # x values of beams
                    beams[i_element, i_end, 0] = np.sum(bay_widths[:i_beam + i_end])
                i_element = i_element + 1
        key = 'beam_geometry'
        hf_group.create_dataset(key, data=beams)

        # store the joint locations
        joints_x = np.array([np.sum(bay_widths[:i_beam]) for i_beam in range(n_bays + 1)])
        # joints_y = np.flip(np.array([np.sum(story_heights[:i_story+1]) for i_story in range(n_stories)]), axis=0)
        joints_y = np.array([np.sum(story_heights[:i_story + 1]) for i_story in range(n_stories)])
        joints_y = np.insert(joints_y, 0, 0, axis=0)  # add the hinge at column base
        [joints_x, joints_y] = np.meshgrid(joints_x, joints_y)
        key = 'joint_locations/joints_x'
        hf_group.create_dataset(key, data=joints_x)
        key = 'joint_locations/joints_y'
        hf_group.create_dataset(key, data=joints_y)


def store_hinge_rotations(results_filename, building_group, model_path):
    ## Creates arrays with the rotation key points of the hinge backbone curves. Stores it in
    # the group with name hf_group inside the h5 file results_filename
    # INPUTS
    #       results_filename: name of h5 file to store frame information
    #       building_group: name of the group to store building metadata
    #       model_path: path to locate the OpenSees tcl file with the material definitions

    with h5py.File(results_filename, 'r+') as hf:
        # Read frame geometry basics
        building_group = hf[building_group]
        n_stories = building_group.attrs['n_stories']
        n_bays = building_group.attrs['n_bays']

        # Open the model file and save in a list of string
        with open(model_path) as model_file:
            model_str = model_file.readlines()
        model_file.close()

        # Create arrays to store the joint rotation capacities
        rotation_yield_positive = np.zeros([n_stories + 1, n_bays + 1, 4, 1])
        rotation_cap_positive = np.zeros([n_stories + 1, n_bays + 1, 4, 1])
        rotation_ultimate_positive = np.zeros([n_stories + 1, n_bays + 1, 4, 1])

        rotation_yield_negative = np.zeros([n_stories + 1, n_bays + 1, 4, 1])
        rotation_cap_negative = np.zeros([n_stories + 1, n_bays + 1, 4, 1])
        rotation_ultimate_negative = np.zeros([n_stories + 1, n_bays + 1, 4, 1])

        # Look for all 4 hinges around each joint
        for floorIdx in range(1, (n_stories + 2)):
            for colIdx in range(1, (n_bays + 2)):
                j = 0
                for eleSide in [1, 2]:
                    # print(eleSide)
                    # eleSide =   # For beams: 1 left hinge in right/current bay and 2 for right hinge in left bay
                    # For columns: 1 is bottom hinge and 2 is top hinge
                    for eleTypeIdx in [2, 3]:
                        #       print(eleTypeIdx)
                        # eleTypeIdx =  # 2 for beams and 3 for columns

                        # Create hinge element tag
                        bayIdx = colIdx + 1 - eleSide
                        storyIdx = floorIdx - 2 + eleSide
                        if floorIdx == 1 and j == 2:
                            # Base column
                            hingeLabel = 300000 + floorIdx * 1000 + colIdx * 10 + 1
                        elif bayIdx > 0 and bayIdx <= n_bays and eleTypeIdx == 3 and floorIdx > 1:
                            # Beam hinges
                            hingeLabel = 200000 + floorIdx * 1000 + bayIdx * 10 + eleSide
                        elif storyIdx <= n_stories and eleTypeIdx == 2 and floorIdx > 1:
                            # Column hinges
                            hingeLabel = 300000 + storyIdx * 1000 + colIdx * 10 + 1
                        else:
                            hingeLabel = 0

                        if hingeLabel != 0:
                            # Find line with the hinge we are looking for

                            # For original models
                            key = 'CreateIbarraMaterial     ' + str(hingeLabel)

                            # For modified models
                            #key = 'CreateIbarraMaterial ' + str(hingeLabel)

                            res = [i for i in model_str if key in i]
                            matLine = str(res)

                            # Extract hinge rotation capacity
                            matAssignments = matLine.split()

                            EIeff = float(matAssignments[3])
                            myPos = float(matAssignments[4])
                            myNeg = float(matAssignments[5])
                            thetaCapPos = float(matAssignments[7])
                            thetaCapNeg = float(matAssignments[8])
                            stiffFactor1 = 11
                            eleLength = float(matAssignments[15])
                            thetaPC = float(matAssignments[9])

                            thetaCap = (thetaCapPos - thetaCapNeg) / 2

                            elstk = stiffFactor1 * ((6 * EIeff) / eleLength)

                            rotation_yield_positive[floorIdx-1, colIdx-1, j, 0] = myPos / elstk
                            rotation_cap_positive[floorIdx-1, colIdx-1, j, 0] = myPos / elstk + thetaCap
                            rotation_ultimate_positive[floorIdx-1, colIdx-1, j, 0] = myPos / elstk + thetaCap + thetaPC

                            rotation_yield_negative[floorIdx-1, colIdx-1, j, 0] = myNeg / elstk
                            rotation_cap_negative[floorIdx-1, colIdx-1, j, 0] = myNeg / elstk - thetaCap
                            rotation_ultimate_negative[floorIdx-1, colIdx-1, j, 0] = myNeg / elstk - thetaCap - thetaPC

                        j = j + 1

        # Save in h5 file's building_group
        key = 'hinge_yield_rotation_positive'
        building_group.create_dataset(key, data=rotation_yield_positive)
        key = 'hinge_cap_rotation_positive'
        building_group.create_dataset(key, data=rotation_cap_positive)
        key = 'hinge_ultimate_rotation_positive'
        building_group.create_dataset(key, data=rotation_ultimate_positive)

        key = 'hinge_yield_rotation_negative'
        building_group.create_dataset(key, data=rotation_yield_negative)
        key = 'hinge_cap_rotation_negative'
        building_group.create_dataset(key, data=rotation_cap_negative)
        key = 'hinge_ultimate_rotation_negative'
        building_group.create_dataset(key, data=rotation_ultimate_negative)

def store_hinge_rotations_modelUQ(results_filename, building_group, model_path, model_id):
    ## Creates arrays with the rotation key points of the hinge backbone curves. Stores it in
    # the group with name hf_group inside the h5 file results_filename
    # INPUTS
    #       results_filename: name of h5 file to store frame information
    #       building_group: name of the group to store building metadata
    #       model_path: path to locate the OpenSees tcl file with the material definitions
    #       model_id: string with name of the group to save hinge properties

    with h5py.File(results_filename, 'r+') as hf:
        # Read frame geometry basics
        building_group = hf[building_group]
        n_stories = building_group.attrs['n_stories']
        n_bays = building_group.attrs['n_bays']

        # Open the model file and save in a list of string
        with open(model_path) as model_file:
            model_str = model_file.readlines()
        model_file.close()

        # Create arrays to store the joint rotation capacities
        rotation_yield_positive = np.zeros([n_stories + 1, n_bays + 1, 4, 1])
        rotation_cap_positive = np.zeros([n_stories + 1, n_bays + 1, 4, 1])
        rotation_ultimate_positive = np.zeros([n_stories + 1, n_bays + 1, 4, 1])

        rotation_yield_negative = np.zeros([n_stories + 1, n_bays + 1, 4, 1])
        rotation_cap_negative = np.zeros([n_stories + 1, n_bays + 1, 4, 1])
        rotation_ultimate_negative = np.zeros([n_stories + 1, n_bays + 1, 4, 1])

        # Look for all 4 hinges around each joint
        for floorIdx in range(1, (n_stories + 2)):
            for colIdx in range(1, (n_bays + 2)):
                j = 0
                for eleSide in [1, 2]:
                    # print(eleSide)
                    # eleSide =   # For beams: 1 left hinge in right/current bay and 2 for right hinge in left bay
                    # For columns: 1 is bottom hinge and 2 is top hinge
                    for eleTypeIdx in [2, 3]:
                        #       print(eleTypeIdx)
                        # eleTypeIdx =  # 2 for beams and 3 for columns

                        # Create hinge element tag
                        bayIdx = colIdx + 1 - eleSide
                        storyIdx = floorIdx - 2 + eleSide
                        if floorIdx == 1 and j == 2:
                            # Base column
                            hingeLabel = 300000 + floorIdx * 1000 + colIdx * 10 + 1
                        elif bayIdx > 0 and bayIdx <= n_bays and eleTypeIdx == 3 and floorIdx > 1:
                            # Beam hinges
                            hingeLabel = 200000 + floorIdx * 1000 + bayIdx * 10 + eleSide
                        elif storyIdx <= n_stories and eleTypeIdx == 2 and floorIdx > 1:
                            # Column hinges
                            hingeLabel = 300000 + storyIdx * 1000 + colIdx * 10 + 1
                        else:
                            hingeLabel = 0

                        if hingeLabel != 0:
                            # Find line with the hinge we are looking for

                            # For original models
                            # key = 'CreateIbarraMaterial     ' + str(hingeLabel)

                            # For modified models
                            key = 'CreateIbarraMaterial ' + str(hingeLabel)

                            res = [i for i in model_str if key in i]
                            matLine = str(res)

                            # Extract hinge rotation capacity
                            matAssignments = matLine.split()

                            EIeff = float(matAssignments[3])
                            myPos = float(matAssignments[4])
                            myNeg = float(matAssignments[5])
                            thetaCapPos = float(matAssignments[7])
                            thetaCapNeg = float(matAssignments[8])
                            stiffFactor1 = 11
                            eleLength = float(matAssignments[15])
                            thetaPC = float(matAssignments[9])

                            thetaCap = (thetaCapPos - thetaCapNeg) / 2

                            elstk = stiffFactor1 * ((6 * EIeff) / eleLength)

                            rotation_yield_positive[floorIdx-1, colIdx-1, j, 0] = myPos / elstk
                            rotation_cap_positive[floorIdx-1, colIdx-1, j, 0] = myPos / elstk + thetaCap
                            rotation_ultimate_positive[floorIdx-1, colIdx-1, j, 0] = myPos / elstk + thetaCap + thetaPC

                            rotation_yield_negative[floorIdx-1, colIdx-1, j, 0] = myNeg / elstk
                            rotation_cap_negative[floorIdx-1, colIdx-1, j, 0] = myNeg / elstk - thetaCap
                            rotation_ultimate_negative[floorIdx-1, colIdx-1, j, 0] = myNeg / elstk - thetaCap - thetaPC

                        j = j + 1

        # Save in h5 file's building_group
        key = 'model_' + model_id + '/hinge_yield_rotation_positive'
        building_group.create_dataset(key, data=rotation_yield_positive)
        key = 'model_' + model_id + '/hinge_cap_rotation_positive'
        building_group.create_dataset(key, data=rotation_cap_positive)
        key = 'model_' + model_id + '/hinge_ultimate_rotation_positive'
        building_group.create_dataset(key, data=rotation_ultimate_positive)

        key = 'model_' + model_id + '/hinge_yield_rotation_negative'
        building_group.create_dataset(key, data=rotation_yield_negative)
        key = 'model_' + model_id + '/hinge_cap_rotation_negative'
        building_group.create_dataset(key, data=rotation_cap_negative)
        key = 'model_' + model_id + '/hinge_ultimate_rotation_negative'
        building_group.create_dataset(key, data=rotation_ultimate_negative)


def collect_gm_metadata(gm_files, results_filename, gm_group):

    [gm_metadata_file,
     gm_sa_avg_file,
     gm_sa_t1_file,
     gm_duration_file,
     gm_spectra_file,
     gm_acc_folder] = gm_files

    gm_metadata = pd.read_csv(gm_metadata_file, sep='\t')

    with open(gm_sa_t1_file, 'r') as file:
        im = file.read().splitlines()
        im = [float(x.strip()) for x in im]
        gm_metadata['Unscaled Sa(T1)'] = im

    with open(gm_sa_avg_file, 'r') as file:
        im = file.read().splitlines()
        im = [float(x.strip()) for x in im]
        gm_metadata['Unscaled Sa_avg'] = im

    sa_t1 = gm_metadata['Unscaled Sa(T1)']
    sa_avg = gm_metadata['Unscaled Sa_avg']

    sa_ratio = sa_t1 / sa_avg
    gm_metadata['Sa_ratio'] = sa_ratio

    durations = pd.read_csv(gm_duration_file, sep='\t')
    gm_metadata['Duration_5-75'] = durations[' Ds575 '].values

    response_spectra = pd.read_csv(gm_spectra_file, sep='\t', index_col='Record ')
    new_index = [x.strip() for x in response_spectra.index]

    # Removes file ext if need be
    if "." in new_index[0]:
        new_index = [x[0:-4] for x in new_index]

    gm_ids = new_index
    gm_metadata['id'] = gm_ids
    gm_metadata.set_index('id', inplace=True)
    response_spectra['id'] = gm_ids

    new_col = [x.strip() for x in response_spectra.columns]
    response_spectra.set_index('id', inplace=True)
    response_spectra.columns = new_col[0:len(response_spectra.columns)]
    periods = [float(x[3:-2]) for x in new_col[0:len(response_spectra.columns)]]

    key = '/ground_motion_records/gm_response_spectra'
    response_spectra.to_hdf(results_filename, key=key)

    with h5py.File(results_filename, 'r+') as hf:
        gm_group = hf[gm_group]
        gm_group['gm_response_spectra'].attrs['periods'] = periods

        # store each ground motion
        for gm_id in gm_ids:
            gm_record_group = gm_group.create_group(gm_id)

            # gm_record_group.attrs['rsn'] = gm_metadata.loc[gm_id, 'RSN']
            gm_record_group.attrs['event'] = gm_metadata.loc[gm_id, 'eventName']
            gm_record_group.attrs['date'] = gm_metadata.loc[gm_id, 'Date']
            gm_record_group.attrs['station'] = gm_metadata.loc[gm_id, 'Station']
            gm_record_group.attrs['magnitude'] = gm_metadata.loc[gm_id, 'M']
            gm_record_group.attrs['r_rup'] = gm_metadata.loc[gm_id, 'Rup']
            gm_record_group.attrs['r_jb'] = gm_metadata.loc[gm_id, 'Rjb']
            # gm_record_group.attrs['vs30'] = gm_metadata.loc[gm_id, 'Vs30']
            # gm_record_group.attrs['region'] = gm_metadata.loc[gm_id, 'region']
            # gm_record_group.attrs['fault_type'] = gm_metadata.loc[gm_id, 'Fault_Type']
            # gm_record_group.attrs['component'] = gm_metadata.loc[gm_id, 'Component']
            gm_record_group.attrs['unscaled_sa_t1'] = gm_metadata.loc[gm_id, 'Unscaled Sa(T1)']
            gm_record_group.attrs['unscaled_sa_avg'] = gm_metadata.loc[gm_id, 'Unscaled Sa_avg']

            # acceleration time history
            if gm_id[0:2] == 'GM':
                gm_acc_file = posixpath.join(gm_acc_folder, gm_id + '.txt')
            else:
                gm_acc_file = posixpath.join(gm_acc_folder, gm_id + '.AT2')

            with open(gm_acc_file, 'r') as file:
                acc = np.array([float(x) for x in file.read().splitlines()])
            dset = gm_record_group.create_dataset('acceleration_time_history', data=acc)
            dset.attrs['n_pts'] = len(acc)
            dset.attrs['dt'] = gm_metadata.loc[gm_id, 'dt']

            # response spectrum
            spectrum = response_spectra.loc[gm_id]
            dset = gm_record_group.create_dataset('response_spectrum', data=spectrum)
            dset.attrs['periods'] = periods

    return gm_metadata


def collect_ida_results_sf(ida_folder, gm_metadata, results_filename, ida_results_group):
    gm_ids = gm_metadata.index
    n_gms = len(gm_ids)

    collapse_values = np.zeros((n_gms, 3))

    # loop through results for each ground motion
    for i in range(n_gms):
        gm_id = gm_ids[i]
        gm_number = int(gm_id[2:])
        idx_number = gm_number - 1

        # read the ida curve
        ida_intensities_file = posixpath.join(ida_folder, gm_id + './ida_curve.txt')
        ida_curve = pd.read_csv(ida_intensities_file, sep='\t', header=None,
                                names=['Scale Factor', 'Story Drift Ratio (max)'])
        # add the collapse point
        last_intensity = ida_curve.iloc[-1:].values[0][0]
        sa_t1_col = last_intensity + 0.01
        ida_curve.loc[len(ida_curve)] = [sa_t1_col, 0.2]
        # add intensities as Sa_avg and scale factors
        sf = ida_curve['Scale Factor'].values
        sa_t1 = sf * gm_metadata['Unscaled Sa(T1)'][idx_number]
        ida_curve.insert(loc=0, column='Sa(T1)', value=sa_t1)
        sa_avg = sf * gm_metadata['Unscaled Sa_avg'][idx_number]
        ida_curve.insert(loc=2, column='Sa_avg', value=sa_avg)

        # save ida curve
        key = ida_results_group + '/' + gm_id + '/ida_curve'
        ida_curve.to_hdf(results_filename, key=key)

        # store collapse values
        collapse_values[idx_number, :] = [sf[-1], sa_t1[-1], sa_avg[-1]]

    # save collapse intensities
    collapse_intensities = pd.DataFrame(collapse_values, index=gm_ids, columns=['Scale Factor', 'Sa(T1)', 'Sa_avg'])
    key = ida_results_group + '/collapse_intensities'
    collapse_intensities.to_hdf(results_filename, key=key)

    # save collapse fragilities
    collapse_fragilities = pd.DataFrame(columns=['Median', 'Beta'], dtype='float64')
    for im in ['Sa(T1)', 'Sa_avg']:
        collapse_fragilities.loc[im, :] = compute_ida_fragility(collapse_intensities[im], plot=True)
    key = ida_results_group + '/collapse_fragilities'
    collapse_fragilities.to_hdf(results_filename, key=key)


def collect_ida_results(ida_folder, gm_metadata, results_filename, ida_results_group):
    gm_ids = gm_metadata.index
    n_gms = len(gm_ids)

    collapse_values = np.zeros((n_gms, 3))

    # loop through results for each ground motion
    for i in range(n_gms):
        gm_id = gm_ids[i]
        # gm_number = int(gm_id[2:])
        idx_number = i #gm_number - 1

        # read the ida curve
        ida_intensities_file = posixpath.join(ida_folder, gm_id + './ida_curve.txt')
        ida_curve = pd.read_csv(ida_intensities_file, sep='\t', header=None,
                                names=['Sa(T1)', 'Story Drift Ratio (max)'])
        # add the collapse point
        last_intensity = ida_curve.iloc[-1:].values[0][0]
        sa_t1_col = last_intensity + 0.01
        ida_curve.loc[len(ida_curve)] = [sa_t1_col, 0.2]
        # add intensities as Sa_avg and scale factors
        sa_t1 = ida_curve['Sa(T1)'].values
        sf = sa_t1 / gm_metadata['Unscaled Sa(T1)'][idx_number]
        ida_curve.insert(loc=0, column='Scale Factor', value=sf)
        sa_avg = sf * gm_metadata['Unscaled Sa_avg'][idx_number]
        ida_curve.insert(loc=2, column='Sa_avg', value=sa_avg)

        # save ida curve
        key = ida_results_group + '/' + gm_id + '/ida_curve'
        ida_curve.to_hdf(results_filename, key=key)

        # store collapse values
        collapse_values[idx_number, :] = [sf[-1], sa_t1[-1], sa_avg[-1]]

    # save collapse intensities
    collapse_intensities = pd.DataFrame(collapse_values, index=gm_ids, columns=['Scale Factor', 'Sa(T1)', 'Sa_avg'])
    key = ida_results_group + '/collapse_intensities'
    collapse_intensities.to_hdf(results_filename, key=key)

    # save collapse fragilities
    collapse_fragilities = pd.DataFrame(columns=['Median', 'Beta'], dtype='float64')
    for im in ['Sa(T1)', 'Sa_avg']:
        collapse_fragilities.loc[im, :] = compute_ida_fragility(collapse_intensities[im], plot=True)
    key = ida_results_group + '/collapse_fragilities'
    collapse_fragilities.to_hdf(results_filename, key=key)


def collect_ida_results_not_finished(ida_folder, gm_metadata, results_filename, ida_results_group):
    gm_ids = gm_metadata.index
    n_gms = len(gm_ids)

    collapse_values = np.zeros((n_gms, 3))

    # loop through results for each ground motion
    for i in range(n_gms):
        gm_id = gm_ids[i]
        gm_number = int(gm_id[2:])
        idx_number = gm_number - 1

        # read the ida curve from completed scales file
        ida_intensities_file = posixpath.join(ida_folder, gm_id + './tolerance_note.txt')
        ida_curve_temp = pd.read_csv(ida_intensities_file, sep='\t', header=None,
                                names=['Sa(T1)', 'Story Drift Ratio (max)', '', '', 'max SDR before collapse', 'convergence', ''])

        sa_t1 = ida_curve_temp['Sa(T1)'].values
        sdr = ida_curve_temp['Story Drift Ratio (max)'].values
        sdr[sdr == np.inf] = 0.10
        order = np.argsort(sdr)
        sdr = sdr[order]
        sa_t1 = sa_t1[order]

        ida_curve = pd.DataFrame(data=np.vstack((sa_t1, sdr)).reshape(len(sdr), 2), columns=['Sa(T1)', 'Story Drift Ratio (max)'])

        # add intensities as Sa_avg and scale factors
        sf = sa_t1 / gm_metadata['Unscaled Sa(T1)'][idx_number]
        ida_curve.insert(loc=0, column='Scale Factor', value=sf)
        sa_avg = sf * gm_metadata['Unscaled Sa_avg'][idx_number]
        ida_curve.insert(loc=2, column='Sa_avg', value=sa_avg)

        # save ida curve
        key = ida_results_group + '/' + gm_id + '/ida_curve'
        ida_curve.to_hdf(results_filename, key=key)

        # store collapse values
        collapse_values[idx_number, :] = [sf[-1], sa_t1[-1], sa_avg[-1]]

    # save collapse intensities
    collapse_intensities = pd.DataFrame(collapse_values, index=gm_ids, columns=['Scale Factor', 'Sa(T1)', 'Sa_avg'])
    key = ida_results_group + '/collapse_intensities'
    collapse_intensities.to_hdf(results_filename, key=key)

    # save collapse fragilities
    collapse_fragilities = pd.DataFrame(columns=['Median', 'Beta'], dtype='float64')
    for im in ['Sa(T1)', 'Sa_avg']:
        collapse_fragilities.loc[im, :] = compute_ida_fragility(collapse_intensities[im], plot=True)
    key = ida_results_group + '/collapse_fragilities'
    collapse_fragilities.to_hdf(results_filename, key=key)



def collect_ida_results_modelUQ(ida_folder, gm_metadata, results_filename, ida_results_group, model_id):
    gm_ids = gm_metadata.index
    n_gms = len(gm_ids)

    collapse_values = np.zeros((n_gms, 3))

    # loop through results for each ground motion
    for i in range(n_gms):
        gm_id = gm_ids[i]
        gm_number = int(gm_id[2:])
        idx_number = gm_number - 1

        # read the ida curve
        ida_intensities_file = posixpath.join(ida_folder, 'AnalysisResult_' + model_id, 'IDA', 'FEMAP695', gm_id + './ida_curve.txt')
        ida_curve = pd.read_csv(ida_intensities_file, sep='\t', header=None,
                                names=['Sa(T1)', 'Story Drift Ratio (max)'])
        # add the collapse point
        last_intensity = ida_curve.iloc[-1:].values[0][0]
        sa_t1_col = last_intensity + 0.01
        ida_curve.loc[len(ida_curve)] = [sa_t1_col, 0.1]
        # add intensities as Sa_avg and scale factors
        sa_t1 = ida_curve['Sa(T1)'].values
        sf = sa_t1 / gm_metadata['Unscaled Sa(T1)'][idx_number]
        ida_curve.insert(loc=0, column='Scale Factor', value=sf)
        sa_avg = sf * gm_metadata['Unscaled Sa_avg'][idx_number]
        ida_curve.insert(loc=2, column='Sa_avg', value=sa_avg)

        # save ida curve
        key = ida_results_group + '/' + gm_id + '/ida_curve/model_' + model_id
        ida_curve.to_hdf(results_filename, key=key)

        # store collapse values
        collapse_values[idx_number, :] = [sf[-1], sa_t1[-1], sa_avg[-1]]

    # save collapse intensities
    collapse_intensities = pd.DataFrame(collapse_values, index=gm_ids, columns=['Scale Factor', 'Sa(T1)', 'Sa_avg'])
    key = ida_results_group + '/collapse_intensities' + '/model_' + model_id
    collapse_intensities.to_hdf(results_filename, key=key)

    # save collapse fragilities
    collapse_fragilities = pd.DataFrame(columns=['Median', 'Beta'], dtype='float64')
    for im in ['Sa(T1)', 'Sa_avg']:
        collapse_fragilities.loc[im, :] = compute_ida_fragility(collapse_intensities[im], plot=True)
    key = ida_results_group + '/collapse_fragilities' + '/model_' + model_id
    collapse_fragilities.to_hdf(results_filename, key=key)


def collect_ida_polarity(gm_scale_folder, results_filename, ida_results_group):
    folder_path = gm_scale_folder.split('AnalysisResult/IDA')
    polarity_folder = folder_path[0] + 'GroundMotion' + folder_path[1]
    polarity_file = posixpath.join(polarity_folder, 'GMafterInfo.txt')

    polarity = pd.read_csv(polarity_file, sep=' ', header=None,
                           names=['Record Number', 'Record Filename', 'dt', 'Polarity', 'empty'])
    polarity['id'] = ['GM' + str(i) for i in polarity['Record Number']]
    polarity.set_index('id', inplace=True)
    polarity = polarity[['Polarity']].copy()

    key = ida_results_group + '/aftershock_polarity'
    polarity.to_hdf(results_filename, key=key)


def collect_ida_time_history(ida_folder, gm_metadata, results_filename, ida_results_group):
    gm_ids = gm_metadata.index
    n_gms = len(gm_ids)

    # loop through results for each ground motion
    for i in range(n_gms):
        gm_id = gm_ids[i]

        key = ida_results_group + '/' + gm_id + '/ida_curve'
        ida_curve = pd.read_hdf(results_filename, key=key)
        im_list = ida_curve['Sa(T1)'].values
        collapse_intensity = im_list[-1]

        aftershock_gm_group = ida_results_group + '/' + gm_id

        with h5py.File(results_filename, 'r+') as hf:
            aftershock_gm_group = hf[aftershock_gm_group]
            time_history_group = aftershock_gm_group.create_group('time_histories')

            building_group = hf['/building_metadata']
            n_stories = building_group.attrs['n_stories']

            for x in [0.25, 0.5, 0.75, 1.0]:
                im = x * collapse_intensity
                idx = np.searchsorted(im_list, im)
                im = '{0:.2f}'.format(im_list[idx])
                im_name = 'Sa(T1)=' + im

                if not im_name in list(time_history_group.keys()):
                    aftershock_gm_folder = ida_folder + '/' + gm_id + './Scale_' + im
                    th_intensity_group = time_history_group.create_group(im_name)

                    file_tag = '_disp.out'
                    filename = posixpath.join(aftershock_gm_folder, 'story1' + file_tag)
                    time_series = np.squeeze(pd.read_csv(filename, sep=' ', header=None).iloc[:, 0].values)
                    time_series = time_series[~np.isnan(time_series)]
                    n_pts = len(time_series)

                    th_intensity_group.create_dataset('time_series', data=time_series)

                    edp = 'displacement'
                    edp_name = 'story_displacement_time_history'
                    n_levels = n_stories
                    file_tag = '_disp.out'

                    edp_results = np.zeros((n_levels, n_pts))
                    for i in range(n_levels):
                        filename = posixpath.join(aftershock_gm_folder, 'story' + str(i + 1) + file_tag)
                        time_history = pd.read_csv(filename, sep=' ', header=None)
                        try:
                            edp_results[i, :] = np.squeeze(time_history.iloc[:, -1])
                        except:
                            edp_results[i, :] = np.squeeze(time_history.iloc[:-1, -1])
                    dset = th_intensity_group.create_dataset(edp_name, data=edp_results)
                    dset.attrs['units'] = 'inches'


def compute_ida_fragility(collapse_ims, plot):
    n_gms = len(collapse_ims)

    median = np.exp((1 / n_gms) * np.sum(np.array([np.log(im) for im in collapse_ims])))
    beta = np.sqrt((1 / (n_gms - 1)) * np.sum(np.array([np.square(np.log(im / median)) for im in collapse_ims])))

    if plot:
        plot_ida_fragility(median, beta, collapse_ims)

    return median, beta


def plot_ida_fragility(median, beta, collapse_ims):
    n_gms = len(collapse_ims)

    plt.scatter(np.sort(collapse_ims), np.linspace(100 / n_gms, 100, num=n_gms, endpoint=True), facecolor='none', edgecolor='lightgray')

    y = np.linspace(0.001, 1, num=100)
    x = stats.lognorm(beta, scale=median).ppf(y)

    label = '$IM_{0.5}=$' + '{0:.2f}'.format(median) + ', $\sigma_{ln}=$' + '{0:.2f}'.format(beta)
    plt.plot(x, 100 * y, label=label, zorder=10, linewidth=3)
    plt.legend()
    plt.show()


def compute_truncated_ida_fragility(collapse_ims, truncation_im, plot):

    parameters_0 = compute_ida_fragility(collapse_ims, plot=False)

    [median, beta] = optimize.minimize(ida_log_likelihood, parameters_0, args=(collapse_ims, truncation_im),
                                       method='Nelder-Mead').x

    if plot:
        plot_truncated_ida_fragility(median, beta, collapse_ims, truncation_im)

    return median, beta


def ida_log_likelihood(parameters, collapse_ims, truncation_im):
    [median, beta] = parameters

    n_gms = len(collapse_ims)

    truncated_collapse_ims = collapse_ims[collapse_ims <= truncation_im]

    im_max = np.max(truncated_collapse_ims)
    n_uncollapsed = n_gms - len(truncated_collapse_ims)

    likelihood_uncollapsed = n_uncollapsed * np.log(1 - stats.lognorm(beta, scale=median).cdf(im_max))

    likelihood_collapsed = np.sum([np.log(stats.lognorm(beta, scale=median).pdf(im)) for im in truncated_collapse_ims])

    log_likelihood = - (likelihood_collapsed + likelihood_uncollapsed)

    return log_likelihood


def plot_truncated_ida_fragility(median, beta, collapse_ims, truncation_im):
    n_gms = len(collapse_ims)

    plt.scatter(np.sort(collapse_ims), np.linspace(100 / n_gms, 100, num=n_gms, endpoint=True), facecolor='none',
                edgecolor='lightgray')

    y = np.linspace(0.001, 1, num=100)
    x = stats.lognorm(beta, scale=median).ppf(y)

    label = '$IM_{0.5}=$' + '{0:.2f}'.format(median) + ', $\sigma_{ln}=$' + '{0:.2f}'.format(beta)
    plt.plot(x, 100 * y, label=label, zorder=10, linewidth=3)

    label = '$IM_{truncated}=$' + '{0:.2f}'.format(truncation_im)
    plt.plot([truncation_im] * 2, [0, 100], color='k', label=label)

    plt.legend()
    plt.show()


def collect_msa_results(msa_folder, gm_metadata, results_filename, msa_results_group):
    gm_ids = gm_metadata.index
    n_gms = len(gm_ids)

    stripe_folders = [i for i in os.listdir(msa_folder) if 'STR' in i]
    n_stripes = len(stripe_folders)
    stripe_values = [float(x[3:]) for x in stripe_folders]

    peak_idr_matrix = np.zeros((n_gms, n_stripes))

    # collect collapse matrix for every ground motion in every stripe
    for i in range(n_stripes):
        for j in range(n_gms):
            msa_idr_file = posixpath.join(msa_folder, stripe_folders[i], gm_ids[j] + '/MSA.txt')
            with open(msa_idr_file, 'r') as file:
                peak_idr_matrix[j, i] = float(file.read())

    # save peak idr matrix
    peak_idr_matrix = pd.DataFrame(peak_idr_matrix, index=gm_ids, columns=stripe_values)
    key = msa_results_group + '/peak_idr_matrix'
    peak_idr_matrix.to_hdf(results_filename, key=key)

    # save collapse matrix
    collapse_matrix = peak_idr_matrix >= 0.1
    collapse_matrix = pd.DataFrame(collapse_matrix, index=gm_ids, columns=stripe_values)
    key = msa_results_group + '/collapse_matrix'
    collapse_matrix.to_hdf(results_filename, key=key)

    # save collapses fragilities
    collapse_fragilities = pd.DataFrame(columns=['Median', 'Beta'], dtype='float64')
    collapse_fragilities.loc[0] = compute_msa_fragility(collapse_matrix, plot=True)
    key = msa_results_group + '/collapse_fragilities'
    collapse_fragilities.to_hdf(results_filename, key=key)


def compute_msa_fragility(collapse_matrix, plot):
    stripe_values = collapse_matrix.columns
    collapse_matrix = collapse_matrix.to_numpy()
    [n_gms, _] = collapse_matrix.shape

    # set the initial median
    p_stripes = np.sum(collapse_matrix, axis=0) / n_gms
    p_target = 0.5
    # linear interpolation for the im resulting in p_target collapses
    if np.any(p_stripes >= p_target):
        median_0 = np.interp(p_target, p_stripes, stripe_values)
    # or take the max im value
    else:
        median_0 = stripe_values[-1]

    [median, beta] = optimize.minimize(msa_log_likelihood, [median_0, 0.3], args=(collapse_matrix, stripe_values),
                                       method='Nelder-Mead').x

    if plot:
        [n_gms, _] = collapse_matrix.shape
        percent_collapses = 100 * np.sum(collapse_matrix, axis=0) / n_gms
        plot_msa_fragility(median, beta, stripe_values, percent_collapses)

    return median, beta


def compute_msa_fragility_plist(p_stripes, stripe_values, plot):
    # compute_msa_fragility takes the fraction of collapse cases for a list of stripes and
    # fits a lognormal probability function
    #
    # INPUTS
    #   p_stripes       = list with the fraction of collapses per stri[e
    #   stripe_values   = list of IM value per stripe
    #   plot            = true/false to plot fragility
    #
    # OUTPUTS
    #   median  = median IM of collapse
    #   beta    = log standard deviation of collapse fragility
    #

    # set the initial median
    p_target = 0.5
    # linear interpolation for the im resulting in p_target collapses
    if np.any(p_stripes >= p_target):
        median_0 = np.interp(p_target, p_stripes, stripe_values)
    # or take the max im value
    else:
        median_0 = stripe_values[-1]

    # standard deviation starts with 0.2
    sigma0 = 0.2

    # convergence flag for optimization
    conv_flag = 0

    while not conv_flag and sigma0 < 1.0:
        # Fit fragility using maximimun likelihood
        x0 = [median_0, sigma0]
        res = optimize.minimize(msa_log_likelihood_new, x0, args=(stripe_values, p_stripes),
                                method='Nelder-Mead', options={'maxiter': 10000})
        conv_flag = res.success
        median, beta = res.x
        sigma0 = sigma0 + 0.05

    if plot:
        y = np.linspace(0.001, 1, num=100)
        x = stats.lognorm(beta, scale=median).ppf(y)
        _ = plt.plot(x, y)
        _ = plt.scatter(stripe_values, p_stripes)

    return median, beta


def msa_log_likelihood_new(parameters, stripe_values, p_stripes):
    # msa_log_likelihood computes the maximum likelihood for a lognormal distribution

    [median, beta] = parameters

    # big sampling number
    bignum = 1000

    num_yy = np.around(bignum * p_stripes).reshape((-1, 1))
    n_stripes = len(stripe_values)

    p_stripes = [stats.lognorm(beta, scale=median).cdf(im) for im in stripe_values]
    stripe_likelihoods = np.array([stats.binom(bignum, p_stripes[i]).pmf(num_yy[i]) for i in range(n_stripes)])

    log_likelihood = - np.sum(np.log(stripe_likelihoods))

    return log_likelihood


def msa_log_likelihood(parameters, collapse_matrix, stripe_values):
    [median, beta] = parameters

    n_collapses = np.sum(collapse_matrix, axis=0)
    [n_gms, n_stripes] = collapse_matrix.shape

    p_stripes = [stats.lognorm(beta, scale=median).cdf(im) for im in stripe_values]

    stripe_likelihoods = np.array([stats.binom(n_gms, p_stripes[i]).pmf(n_collapses[i]) for i in range(n_stripes)])

    log_likelihood = - np.sum(np.log(stripe_likelihoods))

    return log_likelihood


def plot_msa_fragility(median, beta, stripe_values, percent_collapses):
    plt.scatter(stripe_values, percent_collapses, facecolor='none', edgecolor='lightgray')

    y = np.linspace(0.001, 1, num=100)
    x = stats.lognorm(beta, scale=median).ppf(y)

    label = '$IM_{0.5}=$' + '{0:.2f}'.format(median) + ', $\sigma_{ln}=$' + '{0:.2f}'.format(beta)
    plt.plot(x, 100 * y, label=label, zorder=10, linewidth=3)
    plt.legend()
    plt.show()


def collect_damaged_results(damaged_folder, gm_metadata, results_filename, damaged_group, building_group, result_type,
                            main_shock_name='', main_shock_sf=[]):
    gm_ids = gm_metadata.index

    if result_type != 'mainshock_edp':

        all_mainshocks = os.listdir(damaged_folder)

        if main_shock_name == '':
            main_shock_name = gm_ids  # ground motions numbered as GM1, GM2, etc
        else:
            folder_idx = 0

        for gm_id in main_shock_name:

            # Retrieve folder names of back to back IDA runs
            gm_id_mainshocks = [i for i in all_mainshocks if gm_id + '_' in i]

            # Recover the scale for the mainshock from the name of the folder
            if main_shock_name == '':
                scales = np.sort([float(x[-6:-3]) for x in gm_id_mainshocks])
            else:
                scales = main_shock_sf

            for scale in scales:
                # Creates group in the H5 file
                if main_shock_name == '':
                    # mainshock scaled based on collapse intensity
                    gm_scale_group = create_damaged_gm_scale_group(results_filename, damaged_group, gm_id, scale,
                                                                   gm_metadata)
                    scale_name = str(scale) + 'Col'
                    gm_scale_folder = posixpath.join(damaged_folder, gm_id + '_' + scale_name)
                else:
                    # mainshock scaled based on SF
                    gm_id_mainshocks_curr = [i for i in gm_id_mainshocks if gm_id + '_SF' + str(scale) in i]
                    print(gm_id_mainshocks_curr)
                    gm_scale_group = create_damaged_gm_scale_group(results_filename, damaged_group, gm_id, scale, '')
                    # gm_scale_folder = posixpath.join(damaged_folder, gm_id_mainshocks[folder_idx])
                    gm_scale_folder = posixpath.join(damaged_folder, gm_id_mainshocks_curr[0])
                    folder_idx += 1

                print(gm_scale_folder)
                if result_type == 'msa_sa_avg':
                    with h5py.File(results_filename, 'r+') as hf:
                        gm_scale_group = hf[gm_scale_group]
                        msa_results_group = gm_scale_group.create_group('msa_sa_avg').name
                    collect_msa_results(gm_scale_folder, gm_metadata, results_filename, msa_results_group)
                elif result_type == 'ida':
                    with h5py.File(results_filename, 'r+') as hf:
                        gm_scale_group = hf[gm_scale_group]
                        ida_results_group = gm_scale_group.create_group('ida').name
                    collect_ida_results(gm_scale_folder, gm_metadata, results_filename, ida_results_group)
                    collect_ida_polarity(gm_scale_folder, results_filename, ida_results_group)
                    # collect_ida_time_history(gm_scale_folder, gm_metadata, results_filename, ida_results_group)
                else:
                    raise ValueError('Add code for result_type.')

    else:

        scales = np.sort([float(x[3:]) for x in os.listdir(damaged_folder)])

        remove_idx = scales == 1.0
        scales = scales[~remove_idx]

        peak_drift_max = pd.DataFrame(index=gm_ids, columns=scales, dtype='float64')
        residual_drift_max = pd.DataFrame(index=gm_ids, columns=scales, dtype='float64')

        for scale in scales:
            scale_name = 'STR' + str(scale) + '0'

            for gm_id in gm_ids:
                gm_scale_group = create_damaged_gm_scale_group(results_filename, damaged_group, gm_id, scale,
                                                               gm_metadata)
                edp_folder = damaged_folder + '/' + scale_name + '/' + gm_id

                with h5py.File(results_filename, 'r+') as hf:
                    gm_scale_group = hf[gm_scale_group]

                    edp_results_group = gm_scale_group.create_group('mainshock_edp')

                    print('Collecting EDPs for ' + str(scale) + 'Col ' + gm_id)
                    collect_mainshock_edp_results(edp_folder, hf, building_group, edp_results_group)
                    peak_drift_max.loc[gm_id, scale] = edp_results_group['peak_story_drift'].attrs[
                        'max_peak_story_drift']
                    residual_drift_max.loc[gm_id, scale] = edp_results_group['residual_story_drift'].attrs[
                        'max_residual_story_drift']

        key = damaged_group + '/peak_story_drift_max'
        peak_drift_max.to_hdf(results_filename, key=key)
        key = damaged_group + '/residual_drift_max'
        residual_drift_max.to_hdf(results_filename, key=key)


def create_damaged_gm_scale_group(results_filename, damaged_group, gm_id, scale, gm_metadata):
    with h5py.File(results_filename, 'r+') as hf:
        damaged_group = hf[damaged_group]

        if gm_id in damaged_group.keys():
            gm_damaged_group = damaged_group[gm_id]
        else:
            gm_damaged_group = damaged_group.create_group(gm_id)

        if gm_metadata != '':
            scale_name = str(scale) + 'Col'
        else:
            scale_name = 'SF_' + str(scale)

        if scale_name in gm_damaged_group.keys():
            gm_scale_group = gm_damaged_group[scale_name]
        else:
            gm_scale_group = gm_damaged_group.create_group(scale_name)

            if gm_metadata != '':
                collapse_intensity = gm_metadata.loc[gm_id, ['Intact Collapse Scale Factor',
                                                             'Intact Collapse Sa(T1)',
                                                             'Intact Collapse Sa_avg']].values

                gm_scale_group.attrs['scale_factor'] = scale * collapse_intensity[0]
                gm_scale_group.attrs['sa_t1'] = scale * collapse_intensity[1]
                gm_scale_group.attrs['sa_avg'] = scale * collapse_intensity[2]

        return gm_scale_group.name


def collect_mainshock_data(gm_folder_path, mainshock_name, results_filename, mainshock_group):
    gm_acc_file = posixpath.join(gm_folder_path, mainshock_name)
    gm_info_file = posixpath.join(gm_folder_path, 'GMmainInfo.txt')

    with h5py.File(results_filename, 'r+') as hf:
        mainshock_group = hf[mainshock_group]
        # store mainshock ground motion

        if mainshock_name in mainshock_group.keys():
            gm_mainshock_group = mainshock_group[mainshock_name]
        else:
            gm_mainshock_group = mainshock_group.create_group(mainshock_name)

        # acceleration time history
        with open(gm_acc_file, 'r') as file:
            acc = np.array([float(x) for x in file.read().splitlines()])
        dset = gm_mainshock_group.create_dataset('acceleration_time_history', data=acc)
        dset.attrs['n_pts'] = len(acc)

        file = open(gm_info_file, 'r')
        content = file.read()
        content_list = content.split("\t")
        file.close
        dset.attrs['dt'] = float(content_list[2])


def collect_mainshock_data(gm_folder_path, mainshock_name, results_filename, mainshock_group):
    gm_acc_file = posixpath.join(gm_folder_path, mainshock_name)
    gm_info_file = posixpath.join(gm_folder_path, 'GMmainInfo.txt')

    with h5py.File(results_filename, 'r+') as hf:
        mainshock_group = hf[mainshock_group]
        # store mainshock ground motion

        if mainshock_name in mainshock_group.keys():
            gm_mainshock_group = mainshock_group[mainshock_name]
        else:
            gm_mainshock_group = mainshock_group.create_group(mainshock_name)

        # acceleration time history
        with open(gm_acc_file, 'r') as file:
            acc = np.array([float(x) for x in file.read().splitlines()])
        dset = gm_mainshock_group.create_dataset('acceleration_time_history', data=acc)
        dset.attrs['n_pts'] = len(acc)

        file = open(gm_info_file, 'r')
        content = file.read()
        content_list = content.split("\t")
        file.close
        dset.attrs['dt'] = float(content_list[2])


def collect_mainshock_edp_results(edp_folder, hf, building_group, edp_results_group):
    #edp_list = ['drift', 'displacement', 'acceleration']
    edp_list = ['drift', 'displacement', 'jointRotations']

    file_tag = '_disp.out'
    filename = posixpath.join(edp_folder, 'story1' + file_tag)
    time_series = np.squeeze(pd.read_csv(filename, sep=' ', header=None).iloc[:, 0].values)
    time_series = time_series[~np.isnan(time_series)]
    n_pts = len(time_series)

    residual_duration = 2.5
    residual_start_time = time_series[-1] - residual_duration
    residual_start_idx = np.searchsorted(time_series, residual_start_time)

    edp_results_group.create_dataset('time_series', data=time_series)

    building_group = hf[building_group]
    n_stories = building_group.attrs['n_stories']	
    n_bays = building_group.attrs['n_bays']

    for edp in edp_list:

        if edp == 'pfa':
            edp_name = 'peak_floor_acceleration'
            n_levels = n_stories  # ground floor does not have a recorder
            file_tag = '_acc_env.out'

            edp_results = np.zeros(n_levels)
            for i in range(n_levels):
                filename = posixpath.join(edp_folder, 'story' + str(i + 1) + file_tag)
                edp_results[i] = pd.read_csv(filename, sep='\t', header=None).iloc[-1]
            dset = edp_results_group.create_dataset(edp_name, data=edp_results)
            dset.attrs['units'] = 'in/s^2'

        elif edp == 'drift':
            edp_name = 'story_drift_time_history'
            n_levels = n_stories
            file_tag = '_drift.out'

            edp_results = np.zeros((n_levels, n_pts))
            for i in range(n_levels):
                filename = posixpath.join(edp_folder, 'story' + str(i + 1) + file_tag)
                time_history = pd.read_csv(filename, sep=' ', header=None)
                try:
                    edp_results[i,:] = np.squeeze(time_history.iloc[:,-1])
                except:
                    edp_results[i, :] = np.squeeze(time_history.iloc[:-1, -1])
            # dset = edp_results_group.create_dataset(edp_name, data=edp_results)

            edp_name = 'peak_story_drift'
            peak_results = np.max(np.abs(edp_results), axis=1)
            dset = edp_results_group.create_dataset(edp_name, data=peak_results)
            dset.attrs['max_peak_story_drift'] = np.max(peak_results)

            edp_name = 'residual_story_drift'
            residual_results = np.median(edp_results[:, residual_start_idx:], axis=1)
            dset = edp_results_group.create_dataset(edp_name, data=residual_results)
            dset.attrs['max_residual_story_drift'] = np.max(np.abs(residual_results))

            # edp_name = 'between_story_drift_time_history'
            # edp_results = np.append(np.zeros((1, n_pts)), edp_results, axis=0)
            # edp_results = edp_results[1:, :] - edp_results[:-1, :]
            # dset = edp_results_group.create_dataset(edp_name, data=edp_results)

            edp_name = 'residual_between_story_drift'
            residual_results = np.median(edp_results[:, residual_start_idx:], axis=1)
            dset = edp_results_group.create_dataset(edp_name, data=residual_results)


        elif edp == 'displacement':
            edp_name = 'story_displacement_time_history'
            n_levels = n_stories
            file_tag = '_disp.out'

            edp_results = np.zeros((n_levels, n_pts))
            for i in range(n_levels):
                filename = posixpath.join(edp_folder, 'story' + str(i + 1) + file_tag)
                time_history = pd.read_csv(filename, sep=' ', header=None)
                try:
                    edp_results[i,:] = np.squeeze(time_history.iloc[:,-1])
                except:
                    edp_results[i, :] = np.squeeze(time_history.iloc[:-1, -1])
            # dset = edp_results_group.create_dataset(edp_name, data=edp_results)
            # dset.attrs['units'] = 'inches'

            edp_name = 'peak_displacement'
            peak_results = np.max(np.abs(edp_results), axis=1)
            dset = edp_results_group.create_dataset(edp_name, data=peak_results)
            dset.attrs['max_peak_displacement'] = np.max(peak_results)
            dset.attrs['units'] = 'inches'

            edp_name = 'residual_displacement'
            residual_results = np.median(edp_results[:, residual_start_idx:], axis=1)
            dset = edp_results_group.create_dataset(edp_name, data=residual_results)
            dset.attrs['max_residual_displacement'] = np.max(np.abs(residual_results))
            dset.attrs['units'] = 'inches'

        elif edp == 'jointRotations':
            # Read and save every joint rotation demand
            # Find length of time series data and initialize results array
            filename = posixpath.join(edp_folder, 'jointRotations40201.out')
            time_series = np.squeeze(pd.read_csv(filename, sep=' ', header=None).iloc[:, 0].values)
            time_series = time_series[~np.isnan(time_series)]
            n_pts = len(time_series)

            hingeRotDemandsTH = np.zeros([n_stories+1,n_bays+1,4,n_pts])
            hingeRotDemandsPeakPos = np.zeros([n_stories+1,n_bays+1,4,1])
            hingeRotDemandsPeakNeg = np.zeros([n_stories + 1, n_bays + 1, 4, 1])
            hingeRotDemandsResidual = np.zeros([n_stories+1,n_bays+1,4,1])

            # Read rotation of hinges at each joint and store in results array
            for floorIdx in range(1, (n_stories + 2)):
                for colIdx in range(1, (n_bays + 2)):

                    if floorIdx >= 2:
                        ## For the air floors
                        # build joint filename to read
                        jointLabel = 40000 + floorIdx*100 + colIdx
                        jointFile = 'jointRotations' + str(jointLabel) + '.out'

                        # read joint file
                        file_path = posixpath.join(edp_folder, jointFile)
                        results = pd.read_csv(file_path, sep=' ', header=None,
                                              names=['bottom', 'right', 'top', 'left', 'center'])

                        # Save as np array in 4 dimensions (floorIdx, colIdx, hingeLocation, time)
                        hingeRotDemandsTH[floorIdx-1,colIdx-1,0,:] = results.bottom
                        hingeRotDemandsTH[floorIdx-1,colIdx-1,1,:] = results.right
                        hingeRotDemandsTH[floorIdx-1,colIdx-1,2,:] = results.top
                        hingeRotDemandsTH[floorIdx-1,colIdx-1,3,:] = results.left

                        hingeRotDemandsPeakPos[floorIdx-1,colIdx-1,0,0] = max(results.bottom)
                        hingeRotDemandsPeakPos[floorIdx-1,colIdx-1,1,0] = max(results.right)
                        hingeRotDemandsPeakPos[floorIdx-1,colIdx-1,2,0] = max(results.top)
                        hingeRotDemandsPeakPos[floorIdx-1,colIdx-1,3,0] = max(results.left)

                        hingeRotDemandsPeakNeg[floorIdx-1,colIdx-1,0,0] = min(results.bottom)
                        hingeRotDemandsPeakNeg[floorIdx-1,colIdx-1,1,0] = min(results.right)
                        hingeRotDemandsPeakNeg[floorIdx-1,colIdx-1,2,0] = min(results.top)
                        hingeRotDemandsPeakNeg[floorIdx-1,colIdx-1,3,0] = min(results.left)

                        hingeRotDemandsResidual[floorIdx-1,colIdx-1,0,0] = np.median(results.bottom[residual_start_idx:])
                        hingeRotDemandsResidual[floorIdx-1,colIdx-1,1,0] = np.median(results.right[residual_start_idx:])
                        hingeRotDemandsResidual[floorIdx-1,colIdx-1,2,0] = np.median(results.top[residual_start_idx:])
                        hingeRotDemandsResidual[floorIdx-1,colIdx-1,3,0] = np.median(results.left[residual_start_idx:])

                    else:
                        ## For the base of the columns
                        # build column base filename to read
                        elementLabel = 6000 + colIdx*10 + 2
                        elementFile = 'columnBase' + str(elementLabel) + '.out'

                        # read joint file
                        file_path = posixpath.join(edp_folder, elementFile)
                        results = pd.read_csv(file_path, sep=' ', header=None, names=['top'])

                        # Save as np array in 4 dimensions (floorIdx, colIdx, hingeLocation, time)
                        hingeRotDemandsTH[floorIdx-1,colIdx-1,2,:] = results.top
                        hingeRotDemandsPeakPos[floorIdx-1,colIdx-1,2,0] = max(results.top)
                        hingeRotDemandsPeakNeg[floorIdx - 1, colIdx - 1, 2, 0] = min(results.top)
                        hingeRotDemandsResidual[floorIdx-1,colIdx-1,2,0] = np.median(results.top[residual_start_idx:])

            # Defines pointer for storing the joint demands
            # edp_name = 'joint_rotations_time_history'
            # dset = edp_results_group.create_dataset(edp_name, data=hingeRotDemandsTH)
            # dset.attrs['units'] = 'rad'
            # dset.attrs['Hinge locations'] = '0: Bottom; 1: Right; 2: Top; 3: Left'

            edp_name = 'peak_joint_rotations_pos'
            dset = edp_results_group.create_dataset(edp_name, data=hingeRotDemandsPeakPos)
            dset.attrs['units'] = 'rad'
            dset.attrs['Hinge locations'] = '0: Bottom; 1: Right; 2: Top; 3: Left'

            edp_name = 'peak_joint_rotations_neg'
            dset = edp_results_group.create_dataset(edp_name, data=hingeRotDemandsPeakNeg)
            dset.attrs['units'] = 'rad'
            dset.attrs['Hinge locations'] = '0: Bottom; 1: Right; 2: Top; 3: Left'

            edp_name = 'residual_joint_rotations'
            dset = edp_results_group.create_dataset(edp_name, data=hingeRotDemandsResidual)
            dset.attrs['units'] = 'rad'
            dset.attrs['Hinge locations'] = '0: Bottom; 1: Right; 2: Top; 3: Left'

        else:
            raise ValueError('define edp results collection method')

####################################################################################################

def get_EDPstory_response(results_folder, n_stories, file, minrdrift=5e-4):
    # INPUTS
    #    results_folder = path to folder with the results of NLRHA
    #    n_stories      = int with the number of stories
    #    file           = list with any group of the following alternatives
    #
    #                     'drift_env: absolute drift envelope 2D np.array (n_story, 1) with scalar per story
    #                     'drift_env_p: positive dirft envelope 2D np.array (n_story, 1) with scalar per story
    #                     'drift_env_n: negative dirft envelope 2D np.array (n_story, 1) with scalar per story
    #
    #                     'acc_env': 2D np.array (n_story, 1) with scalar per story
    #                     'acc_env_p': positive dirft envelope 2D np.array (n_story, 1) with scalar per story
    #                     'acc_env_n': negative dirft envelope 2D np.array (n_story, 1) with scalar per story
    #
    #                     'drift_max: absolute drift envelope (from drift.out files) 2D np.array (n_story, 1) with scalar per story
    #                     'drift_max_p: absolute drift envelope (from drift.out files) 2D np.array (n_story, 1) with scalar per story
    #                     'drift_max_n: absolute drift envelope (from drift.out files) 2D np.array (n_story, 1) with scalar per story
    #
    #                     'rdrift_all': 2D np.array (n_story, 1) with scalar per story keeping the sign
    #                     'rdrift_all_abs: 2D np.array (n_story, 1) with scalar per story
    #                     'rdrift_max': scalar with maximum absolute residual for the building
    #
    #   minrdrift       = float as minimum value of residual drift to consider
    # OUTPUTS
    #    response = np.array with the response desired per story/floor
    #

    # Read results as 1d array
    if 'rdrift' in file or 'drift_max' in file:
        file_save = file
        file = 'drift'
    elif 'drift_env' in file:
        file_save = file
        file = 'drift_env'
    elif 'acc_env' in file:
        file_save = file
        file = 'acc_env'
    else:
        print('File name not supported')

    # Reads the first story/floor
    # (BREAKS DATA COLLECTION IF EMPTY ACC FILES, MEANING THE RHA DID NOT FINISH)
    if file == 'acc_env':
        neg_line = 0
        pos_line = 0
        last_line = 0

        try:
            filepath = posixpath.join(results_folder, 'story' + str(0) + '_' + file + '.out')
            with open(filepath) as f:
                for line in f:
                    try:
                        float(line)  # only place if can convert to float (avoid blank spaces in the file)
                        if line == '\n':
                            last_line = last_line
                        else:
                            last_line = line
                    except:
                        last_line = last_line

                    aux = float(last_line)
                    if neg_line == 0:
                        neg_line = aux
                    elif pos_line == 0:
                        pos_line = aux
                    else:
                        last_line = aux
        except:
            response = 0
            return response

        if '_n' in file_save:
            response = np.array(neg_line)
        elif '_p' in file_save:
            response = np.array(pos_line)
        else:
            response = np.array(last_line)  # read last line as list of strings: ignores first entry (time)

    elif file == 'drift_env':
        neg_line = 0
        pos_line = 0
        last_line = 0

        try:
            filepath = posixpath.join(results_folder, 'story' + str(1) + '_' + file + '.out')
            with open(filepath) as f:
                for line in f:
                    line = line.strip()
                    line = line.split('\n')[0]
                    try:
                        aux = float(line)  # only place if can convert to float (avoid blank spaces in the file)
                    except:
                        aux = 0

                    if neg_line == 0:
                        neg_line = aux
                    elif pos_line == 0:
                        pos_line = aux
                    else:
                        last_line = aux
        except:
            response = 0
            return response

        if '_n' in file_save:
            response = np.array(neg_line)
        elif '_p' in file_save:
            response = np.array(pos_line)
        else:
            response = np.array(last_line)  # read last line as list of strings: ignores first entry (time)

    else:
        # Find story with drift.out file available
        for i in range(n_stories):
            filepath = posixpath.join(results_folder, 'story' + str(i + 1) + '_' + file + '.out')
            try:
                #                 response = pd.read_csv(filepath, header=None, sep=' ')
                #                 response = response.values
                with open(filepath) as f:
                    last_line = 0
                    for line in f:
                        prev_line = last_line
                        if file == 'drift':
                            try:
                                aux = line.split()
                                aux = [float(x) for x in
                                       aux]  # convert to list of floats, empty "line" gives an error here
                                # should be two cell vector because the output must include time or should not be a space
                                if len(aux) == 1 or line == '\n':
                                    last_line = last_line
                                else:
                                    last_line = line
                            except:
                                last_line = last_line
                        else:
                            last_line = line
                    if file == 'drift':
                        aux = prev_line.split()  # read second to last line in a time history since sometimes has errors
                    else:
                        aux = last_line.split()  # read last line as list of strings (drift_env)

                    aux = [float(x) for x in aux]  # convert to list of floats
                    response = np.array(np.abs(aux))
                i_story_worked = i + 1
                break
            except:
                 pass
                # print('MISSING: story' + str(i + 1) + '_' + file + '.out')

        # Check if results on every story are not consistent
        # (usually when the building collapses in the first step)
        if 'i_story_worked' in locals():
            filepath = posixpath.join(results_folder, 'story' + str(i_story_worked) + '_' + file + '.out')
        else:
            print('ERROR COLLECTING DATA FOR:')
            print(results_folder)
            print('')
            return 0

        with open(filepath) as f:
            last_line = 0
            drift_max = 0
            drift_max_p = 0
            drift_max_n = 0
            for line in f:
                prev_line = last_line
                try:
                    aux = line.split()
                    aux = [float(x) for x in aux]  # convert to list of floats, empty "line" gives an error here
                    # should be two cell vector because the output must include time or should not be a space
                    if len(aux) == 1 or line == '\n':
                        last_line = last_line
                    else:
                        last_line = line
                except:
                    last_line = last_line

                aux = last_line.split()
                aux = [float(x) for x in aux]
                drift_max = np.max([np.abs(aux[1]), drift_max])
                drift_max_p = np.max([aux[1], drift_max_p])
                drift_max_n = np.min([aux[1], drift_max_n])

            if file_save == 'drift_max':
                res = np.array(drift_max)
            elif file_save == 'drift_max_p':
                res = np.array(drift_max_p)
            elif file_save == 'drift_max_n':
                res = np.array(drift_max_n)
            elif 'rdrift' in file_save:
                aux = prev_line.split()  # read second to last line in a time history since sometimes has errors
                aux = [float(x) for x in aux]  # convert to list of floats
                if '_abs' in file_save or '_max' in file_save:
                    res = np.array(np.abs(aux[1]))  # ignores first entry (time)
                else:
                    res = np.array(aux[1])  # ignores first entry (time)

            if file_save == 'rdrift_all_abs' or file_save == 'rdrift_max':
                res = max(res, minrdrift)
            elif 'rdrift' in file_save:
                if res < 0:
                    res = min(res, -minrdrift)
                else:
                    res = max(res, minrdrift)
            response = res

    # Reads and append the results of the remaining stories/floors
    if file == 'acc_env':
        for i_floor in range(n_stories):
            filepath = posixpath.join(results_folder, 'story' + str(i_floor) + '_' + file + '.out')
            #             res = np.loadtxt(filepath)
            #             res = pd.read_csv(filepath, header=None, sep=' ')
            #             res = res.values
            with open(filepath) as f:
                neg_line = 0
                pos_line = 0
                last_line = 0
                for line in f:
                    try:
                        float(line)  # only place if can convert to float (avoid blank spaces in the file)
                        # should not be a space
                        if line == '\n':
                            last_line = last_line
                        else:
                            last_line = line
                    except:
                        last_line = last_line

                    aux = float(last_line)
                    if neg_line == 0:
                        neg_line = aux
                    elif pos_line == 0:
                        pos_line = aux
                    else:
                        last_line = aux

                if '_n' in file_save:
                    res = np.array(neg_line)
                elif '_p' in file_save:
                    res = np.array(pos_line)
                else:
                    res = np.array(last_line)  # read last line as list of strings: ignores first entry (time)
            res
            response = np.vstack((response, res))
    else:
        for i_story in range(n_stories - 1):
            i_story = i_story + 1
            filepath = posixpath.join(results_folder, 'story' + str(i_story + 1) + '_' + file + '.out')

            if file == 'drift':
                try:
                    # Check if drift.out file exist (some runs are not producing this output for some stories)
                    # res = pd.read_csv(filepath, header=None, sep=' ')
                    # res = res.values
                    with open(filepath) as f:
                        last_line = 0
                        drift_max = 0
                        drift_max_p = 0
                        drift_max_n = 0
                        for line in f:
                            prev_line = last_line
                            try:
                                aux = line.split()
                                aux = [float(x) for x in
                                       aux]  # convert to list of floats, empty "line" gives an error here
                                # should be two cell vector because the output must include time or should not be a space
                                if len(aux) == 1 or line == '\n':
                                    last_line = last_line
                                else:
                                    last_line = line
                            except:
                                last_line = last_line

                            aux = last_line.split()
                            aux = [float(x) for x in aux]
                            drift_max = np.max([np.abs(aux[1]), drift_max])
                            drift_max_p = np.max([aux[1], drift_max_p])
                            drift_max_n = np.min([aux[1], drift_max_n])

                        if file_save == 'drift_max':
                            res = np.array(drift_max)
                        elif file_save == 'drift_max_p':
                            res = np.array(drift_max_p)
                        elif file_save == 'drift_max_n':
                            res = np.array(drift_max_n)
                        elif 'rdrift' in file_save:
                            aux = prev_line.split()  # read second to last line in a time history since sometimes has errors
                            aux = [float(x) for x in aux]  # convert to list of floats
                            if '_abs' in file_save or '_max' in file_save:
                                res = np.array(np.abs(aux[1]))  # ignores first entry (time)
                            else:
                                res = np.array(aux[1])  # ignores first entry (time)

                    i_story_worked = i_story + 1
                except:
                    # If the drift.out file does not exist for this story, take the previous story
                    print('MISSING: story' + str(
                        i_story + 1) + '_' + file + '.out, TAKING previous story that worked: story' + str(
                        i_story_worked))
                    filepath = posixpath.join(results_folder, 'story' + str(i_story_worked) + '_' + file + '.out')
                    # res = pd.read_csv(filepath, header=None, sep=' ')
                    # res = res.values
                    with open(filepath) as f:
                        last_line = 0
                        drift_max = 0
                        drift_max_p = 0
                        drift_max_n = 0
                        for line in f:
                            prev_line = last_line
                            try:
                                aux = line.split()
                                aux = [float(x) for x in
                                       aux]  # convert to list of floats, empty "line" gives an error here
                                # should be two cell vector because the output must include time or should not be a space
                                if len(aux) == 1 or line == '\n':
                                    last_line = last_line
                                else:
                                    last_line = line
                            except:
                                last_line = last_line

                            aux = last_line.split()
                            aux = [float(x) for x in aux]
                            drift_max = np.max([np.abs(aux[1]), drift_max])
                            drift_max_p = np.max([aux[1], drift_max_p])
                            drift_max_n = np.min([aux[1], drift_max_n])

                        if file_save == 'drift_max':
                            res = np.array(drift_max)
                        elif file_save == 'drift_max_p':
                            res = np.array(drift_max_p)
                        elif file_save == 'drift_max_n':
                            res = np.array(drift_max_n)
                        elif 'rdrift' in file_save:
                            aux = prev_line.split()  # read second to last line in a time history since sometimes has errors
                            aux = [float(x) for x in aux]  # convert to list of floats
                            if '_abs' in file_save or '_max' in file_save:
                                res = np.array(np.abs(aux[1]))  # ignores first entry (time)
                            else:
                                res = np.array(aux[1])  # ignores first entry (time)

            else:  # drift_env.out

                # res = pd.read_csv(filepath, header=None, sep=' ')
                # res = res.values
                neg_line = 0
                pos_line = 0
                last_line = 0
                with open(filepath) as f:
                    for line in f:
                        line = line.strip()
                        line = line.split('\n')[0]
                        try:
                            aux = float(line)  # only place if can convert to float (avoid blank spaces in the file)
                        except:
                            aux = 0

                        if neg_line == 0:
                            neg_line = aux
                        elif pos_line == 0:
                            pos_line = aux
                        else:
                            last_line = aux

                    if '_n' in file_save:
                        res = np.array(neg_line)
                    elif '_p' in file_save:
                        res = np.array(pos_line)
                    else:
                        res = np.array(last_line)  # read last line as list of strings: ignores first entry (time)

            # Minimum RID = minrdrift to avoid numerical issues when pelicun fits a probabilistic model
            if file_save == 'rdrift_all_abs' or file_save == 'rdrift_max':
                res = max(res, minrdrift)
            elif 'rdrift' in file_save:
                if res < 0:
                    res = min(res, -minrdrift)
                else:
                    res = max(res, minrdrift)

            # Save the maximum across floors or results per floor
            if file_save == 'rdrift_max':
                response = np.max([response, res])
            else:
                response = np.vstack((response, res))

    return response


def get_DSC(results_folder, n_connections):
    # Read response for each flange in each connection and interpretes the damage state of the connection per FEMAP58
    #
    # INPUTS
    #    results_folder = path to folder with the results of NLRHA
    #
    # OUTPUTS
    #    DSC            = 1D np.array with the DS for each connection starting with all the
    #                     left connection in each beam and then all the right connections.
    #                     DS0: No damage
    #                     DS1: Fracture of the bottom flange
    #                     DS2: Fracture of the top flange
    #                     DS3: Fracture of both flanges
    #

    # Left connection
    file = 'frac_LB'
    filepath = posixpath.join(results_folder, file + '.out')
    with open(filepath) as f:
        last_line = 0
        prev_line = 0
        pp_line = 0
        for line in f:
            ppp_line = pp_line
            pp_line = prev_line
            prev_line = last_line
            if line == '\n':
                last_line = last_line
            else:
                last_line = line
        aux = last_line.split()  # read second to last line in a time history since sometimes has errors
        # only keeps set of results that include all connections (considering last 4 lines)
        if len(aux) != n_connections:
            aux = prev_line.split()  # read second to last line in a time history since sometimes has errors
            if len(aux) != n_connections:
                aux = pp_line.split()  # read second to last line in a time history since sometimes has errors
                if len(aux) != n_connections:
                    aux = ppp_line.split()  # read second to last line in a time history since sometimes has errors
        aux = [float(x) for x in aux]  # convert to list of floats
        frac_LB = np.abs(aux)

    file = 'frac_LT'
    filepath = posixpath.join(results_folder, file + '.out')
    with open(filepath) as f:
        last_line = 0
        prev_line = 0
        pp_line = 0
        for line in f:
            ppp_line = pp_line
            pp_line = prev_line
            prev_line = last_line
            if line == '\n':
                last_line = last_line
            else:
                last_line = line
        aux = last_line.split()  # read second to last line in a time history since sometimes has errors
        # only keeps set of results that include all connections (considering last 4 lines)
        if len(aux) != n_connections:
            aux = prev_line.split()  # read second to last line in a time history since sometimes has errors
            if len(aux) != n_connections:
                aux = pp_line.split()  # read second to last line in a time history since sometimes has errors
                if len(aux) != n_connections:
                    aux = ppp_line.split()  # read second to last line in a time history since sometimes has errors
        aux = [float(x) for x in aux]  # convert to list of floats
        frac_LT = np.abs(aux)
    frac_LT = frac_LT * 2

    frac_L = frac_LB + frac_LT

    # Right connection
    file = 'frac_RB'
    filepath = posixpath.join(results_folder, file + '.out')
    with open(filepath) as f:
        last_line = 0
        prev_line = 0
        pp_line = 0
        for line in f:
            ppp_line = pp_line
            pp_line = prev_line
            prev_line = last_line
            if line == '\n':
                last_line = last_line
            else:
                last_line = line
        aux = last_line.split()  # read second to last line in a time history since sometimes has errors
        # only keeps set of results that include all connections (considering last 4 lines)
        if len(aux) != n_connections:
            aux = prev_line.split()  # read second to last line in a time history since sometimes has errors
            if len(aux) != n_connections:
                aux = pp_line.split()  # read second to last line in a time history since sometimes has errors
                if len(aux) != n_connections:
                    aux = ppp_line.split()  # read second to last line in a time history since sometimes has errors
        aux = [float(x) for x in aux]  # convert to list of floats
        frac_RB = np.abs(aux)
    file = 'frac_RT'
    filepath = posixpath.join(results_folder, file + '.out')
    with open(filepath) as f:
        last_line = 0
        prev_line = 0
        pp_line = 0
        for line in f:
            ppp_line = pp_line
            pp_line = prev_line
            prev_line = last_line
            if line == '\n':
                last_line = last_line
            else:
                last_line = line
        aux = last_line.split()  # read second to last line in a time history since sometimes has errors
        # only keeps set of results that include all connections (considering last 4 lines)
        if len(aux) != n_connections:
            aux = prev_line.split()  # read second to last line in a time history since sometimes has errors
            if len(aux) != n_connections:
                aux = pp_line.split()  # read second to last line in a time history since sometimes has errors
                if len(aux) != n_connections:
                    aux = ppp_line.split()  # read second to last line in a time history since sometimes has errors
        aux = [float(x) for x in aux]  # convert to list of floats
        frac_RT = np.abs(aux)
    frac_RT = frac_RT * 2

    frac_R = frac_RB + frac_RT

    DSC = np.hstack((frac_L, frac_R))

    return DSC


def get_DSsplice(results_folder, splice_frac_strain, n_splices):
    # Read response for each flange in each connection and interpretes the damage state of the connection per FEMAP58
    #
    # INPUTS
    #    results_folder     = path to folder with the results of NLRHA
    #    splice_frac_strain = strain limit to judge that fracture occured in the splice
    #
    # OUTPUTS
    #    DSsplice           = 1D np.array with the DS for each connection starting with all the
    #                         left connection in each beam and then all the right connections.
    #                         DS0: No damage
    #                         DS1: Fracture of splice weld in either side
    #

    # Inputs
    file = 'ss_splice'

    filepath = posixpath.join(results_folder, file + '.out')
    with open(filepath) as f:
        max_res = np.zeros(2*n_splices)
        last_line = 0
        for line in f:
            prev_line = last_line
            try:
                aux = line.split()
                aux = [abs(float(x)) for x in aux]  # convert to list of floats, empty "line" gives an error here
                if len(aux) != 2*n_splices: # only keeps set of results that include all splices (stress and strain)
                    last_line = last_line
                else:
                    last_line = line
                    max_res = np.max(max_res, aux, axis=0)
            except:
                last_line = last_line
        # Take the last value
        # aux = prev_line.split()  # read second to last line in a time history since sometimes has errors
        # aux = [float(x) for x in aux]  # convert to list of floats
        # data_1d = np.abs(aux)

        # Take the max value in the history
        data_1d = max_res
    n_cols = len(data_1d)

    # extract strain
    strain_splices = data_1d[1:n_cols:2]

    # get boolean for fracture if excesive strain
    DSsplice = np.zeros(len(strain_splices))
    DSsplice[strain_splices > splice_frac_strain] = 1

    return DSsplice


def collect_gmset_response(stripe_folder_path, beam_list, fracElement, dir_i, spliceElement, splice_list, column_list,
                           drift_out='abs', rdrift_out='max', minrdrift=5e-4, splice_frac_strain=60 * 2 / 29000):
    # Creates a table per stripe with the peak responses per story/floor
    #
    # INPUTS
    #    stripe_folder_path = path to find the results for each ground motions
    #    beam_list          = 2D np.array with 1 or 0 for the beams that exist
    #    fracElement        = true : collects damage state for each connection
    #                         false: skip collection of connection damage state
    #    dir_i              = Code used to identify direction for EDP file format
    #                         1: denotes X
    #                         2: denotes Y
    #    spliceElement      = boolean to collect or not splice damage states
    #    splice_list        = 2D np.array with 1 or 0 for the stories where splices exist
    #    column_list        = list of 2D np.array indicating which columns exist
    #    drift_out          = Method to output peak drift
    #                         -> 'abs' = multiple columns with the peak absolute value for each floor
    #                         -> 'both' = multiple columns with the peak positive an negative value for each floor
    #    rdrift_out         = Method to output residual drift
    #                         -> 'max' = single column with the maximum residual across floors
    #                         -> 'all_abs' = multiple columns with the residual absolute value for each floor
    #                         -> 'all' = multiple columns with the residual for each floor with its sign
    #    minrdrift          = float as minimum value of residual drift to consider
    #    splice_frac_strain = strain limit to judge that fracture occured in the splice
    #
    # OUTPUT
    #    response_matrix    = pd.DataFrame with all the results (columns) for each ground motion (rows) in this stripe
    #

    n_stories, n_bays = beam_list.shape

    gm_ids = os.listdir(stripe_folder_path)
    n_gms = len(gm_ids)

    # Define column names for dataframe
    column_names = []
    column_names.append('EndCriteria')
    if drift_out == 'abs':
        for i_story in range(n_stories):
            column_names.append('1-PID-' + str(i_story + 1) + '-' + str(dir_i))
    else:
        for i_story in range(n_stories):
            column_names.append('1-PIDp-' + str(i_story + 1) + '-' + str(dir_i))
        for i_story in range(n_stories):
            column_names.append('1-PIDn-' + str(i_story + 1) + '-' + str(dir_i))
    if 'all' in rdrift_out:
        for i_story in range(n_stories):
            column_names.append('1-RID-' + str(i_story + 1) + '-' + str(dir_i))
    else:
        column_names.append('1-RID-1' + '-' + str(dir_i))
    for i_floor in range(n_stories):
        column_names.append('1-PFA-' + str(i_floor) + '-' + str(dir_i))
    column_names.append('1-PFA-' + str(n_stories) + '-' + str(dir_i))

    if fracElement:
        for i_floor in range(n_stories):
            for i_beam in range(n_bays):
                if beam_list[i_floor, i_beam] > 0:
                    column_names.append('1-DSC-' + str(i_floor + 1) + '-' + "{0:0=4d}".format(
                        dir_i * 1000 + i_beam * 10 + 1))  # DSC: Damage State Connection (left side)
        for i_floor in range(n_stories):
            for i_beam in range(n_bays):
                if beam_list[i_floor, i_beam] > 0:
                    column_names.append('1-DSC-' + str(i_floor + 1) + '-' + "{0:0=4d}".format(
                        dir_i * 1000 + i_beam * 10 + 2))  # DSC: Damage State Connection (right side)

    if spliceElement:
        n_stories, n_pier = column_list.shape
        for i_pier in range(n_pier):
            for i_story in range(n_stories):
                if splice_list[i_story, i_pier] > 0:
                    # jump if no splice in this column segment. The splice_list array already accounts for setbacks, atriums or missing columns
                    column_names.append('1-DSS-' + str(i_story + 1) + '-' + "{0:0=3d}".format(
                        dir_i * 100 + i_pier))  # DSC: Damage State Splice

    # print(len(column_names))
    # print(column_names)

    # collect response per story/floor as a list of arrays, each entry in the list is on gm
    response = []
    gm_ids = os.listdir(stripe_folder_path)
    removeGMlist = []
    for j in range(n_gms):
        # print(gm_ids[j])
        results_folder = posixpath.join(stripe_folder_path, gm_ids[j])
        pfa_gm = get_EDPstory_response(results_folder, n_stories, 'acc_env')
        rdrift_gm = get_EDPstory_response(results_folder, n_stories, 'rdrift_' + rdrift_out, minrdrift=minrdrift)
        if drift_out == 'abs':
            pid_gm = get_EDPstory_response(results_folder, n_stories, 'drift_env')
            if type(pid_gm) == int:
                pid_gm = get_EDPstory_response(results_folder, n_stories, 'drift_max')
        else:
            pid_gm_p = get_EDPstory_response(results_folder, n_stories, 'drift_env_p')
            if type(pid_gm_p) == int:
                pid_gm_p = get_EDPstory_response(results_folder, n_stories, 'drift_max_p')
            pid_gm_n = get_EDPstory_response(results_folder, n_stories, 'drift_env_n')
            if type(pid_gm_n) == int:
                pid_gm_n = get_EDPstory_response(results_folder, n_stories, 'drift_max_n')
            pid_gm = np.hstack((pid_gm_p.flatten(), pid_gm_n.flatten()))

        if type(pfa_gm) == int:  # did not finish RHA, so skip the ground motion
            print('Did not finish GM (ACC ZERO): ' + results_folder)
            print()
            removeGMlist.append(j)
        elif type(rdrift_gm) == int:  # did not finish RHA, so skip the ground motion
            print('Did not finish GM (RDRIFT ZERO): ' + results_folder)
            print()
            removeGMlist.append(j)
        elif type(pid_gm) == int:  # did not finish RHA, so skip the ground motion
            print('Did not finish GM (DRIFT ZERO): ' + results_folder)
            print()
            removeGMlist.append(j)
        else:
            # pid_gm = get_EDPstory_response(results_folder, n_stories, 'drift_env')
            # try:
            # rdrift_gm = get_EDPstory_response(results_folder, n_stories, 'rdrift_' + rdrift_out, minrdrift=minrdrift)
            # except:
            #     print('ERROR drift for ' + gm_ids[j])

            if fracElement:
                dsc_gm = get_DSC(results_folder, np.sum(beam_list))
                # print(np.sum(beam_list)*2)
                # print(len(dsc_gm.flatten()))
            if spliceElement:
                dssplice_gm = get_DSsplice(results_folder, splice_frac_strain, np.sum(splice_list))
                # print(np.sum(splice_list))
                # print(len(dssplice_gm.flatten()))

            filepath = posixpath.join(results_folder, 'MSA.txt')
            try:
                maxDrift_wcollapse = np.loadtxt(filepath)
            except:
                with open(filepath) as f:
                    for line in f:
                        maxDrift_wcollapse = line.strip()

            if type(maxDrift_wcollapse) == str:
                if  maxDrift_wcollapse == 'Collapsed':
                    endCriteria = 'MaxDrift'
                else:
                    endCriteria = 'nonCollapse'

            else:

                # Find unknown errors
                try:
                    np.max(np.abs(pid_gm))
                except:
                    print(pid_gm)
                    print('UNKNOWN ERROR: ' + results_folder)

                # Define endCriteria
                if (maxDrift_wcollapse > 0.099) and (np.max(np.abs(pid_gm)) > 0.07):
                    # Maximum drift reached during the analysis
                    endCriteria = 'MaxDrift'
                #                 response_gm = np.hstack((endCriteria, np.zeros(pid_gm.shape).flatten(), np.zeros(rdrift_gm.shape).flatten(), \
                #                                          np.zeros(pfa_gm.shape).flatten(), np.zeros(dsc_gm.shape).flatten()))
                elif (maxDrift_wcollapse > 0.099) and (np.max(np.abs(pid_gm)) < 0.07):
                    # Inconvergence
                    endCriteria = 'Inconvergence'
                #                 response_gm = np.hstack((endCriteria, np.zeros(pid_gm.shape).flatten(), np.zeros(rdrift_gm.shape).flatten(), \
                #                                          np.zeros(pfa_gm.shape).flatten(), np.zeros(dsc_gm.shape).flatten()))
                else:
                    # Non collapse
                    endCriteria = 'nonCollapse'
                #                 response_gm = np.hstack((endCriteria, pid_gm.flatten(), rdrift_gm.flatten(), pfa_gm.flatten(), dsc_gm.flatten()))

            if fracElement and spliceElement:
                response_gm = np.hstack(
                    (endCriteria, pid_gm.flatten(), rdrift_gm.flatten(), pfa_gm.flatten(), dsc_gm.flatten(),
                     dssplice_gm.flatten()))
            elif fracElement:
                response_gm = np.hstack(
                    (endCriteria, pid_gm.flatten(), rdrift_gm.flatten(), pfa_gm.flatten(), dsc_gm.flatten()))
            elif spliceElement:
                response_gm = np.hstack(
                    (endCriteria, pid_gm.flatten(), rdrift_gm.flatten(), pfa_gm.flatten(), dssplice_gm.flatten()))
            else:
                try:
                    pid_gm.flatten()
                except:
                    print('pid_gm')
                    print(results_folder)
                try:
                    rdrift_gm.flatten()
                except:
                    print('pid_gm')
                    print(results_folder)
                try:
                    pfa_gm.flatten()
                except:
                    print('pid_gm')
                    print(results_folder)

                response_gm = np.hstack((endCriteria, pid_gm.flatten(), rdrift_gm.flatten(), pfa_gm.flatten()))

            # print(response_gm)
            # print(len(response_gm))
            # print('DONE: ' + gm_ids[j])
            response.append(response_gm)

    # save peak idr matrix
    gm_ids = np.delete(gm_ids, removeGMlist)  # remove the gm that did not finish RHA
    # print(len(response))
    # print(len(gm_ids))
    response_matrix = pd.DataFrame(response, columns=column_names, index=gm_ids)

    return response_matrix


def collect_XandY_response(model_name_all, stripe_folder_all, save_results_folder_all, msa_folders_all, beam_list_x_all,
                           beam_list_y_all, fracElement, spliceElement_all, splice_list_x_all, splice_list_y_all,
                           column_list_x_all,
                           column_list_y_all, minrdrift, splice_frac_strain, case_i):
    # INPUTS
    # Collects the EDP results considering each stripe of each case as an independent job
    #    model_name_all           = list of str with the case name to collect results from
    #    stripe_folder_all        = list of str with the foldername of the stripe to collect results from
    #    save_results_folder_all  = list of str with path to save results
    #    msa_folders_all          = list of list [path_to_results_X, path_to_results_Y]
    #    beam_list_x_all          = list of 2D np.array indicating which beams exist in the X frame
    #    beam_list_Y_all          = list of 2D np.array indicating which beams exist in the Y frame
    #    fracElement              = true  -> collect connections DS
    #                               false -> does NOT collect connections DS
    #
    #    spliceElement_all        = 1 -> collect splice data
    #                               0 -> do not collect splice data
    #    splice_list_x_all        = list of 2D np.array indicating the pier and story with splice in the X frame
    #    splice_list_y_all        = list of 2D np.array indicating the pier and story with splice in the Y frame
    #    column_list_x_all        = list of 2D np.array indicating which columns exist in the X frame
    #    column_list_y_all        = list of 2D np.array indicating which column exist in the X frame
    #    minrdrift                = minimum residual drift to collect
    #    splice_frac_strain       = strain limit for deciding fracture occured on splices
    #

    # Parse case to execute
    model_name = model_name_all[case_i]
    stripe_folder = stripe_folder_all[case_i]
    save_results_folder = save_results_folder_all[case_i]
    msa_folders = msa_folders_all[case_i]
    beam_list_x = beam_list_x_all[case_i]
    beam_list_y = beam_list_y_all[case_i]

    spliceElement = spliceElement_all[case_i]
    column_list_x = column_list_x_all[case_i]
    column_list_y = column_list_y_all[case_i]
    splice_list_x = splice_list_x_all[case_i]
    splice_list_y = splice_list_y_all[case_i]

    # Define path to save results
    results_filename = posixpath.join(save_results_folder, 'EDP_' + model_name + '_' + stripe_folder + '.csv')

    #### collect results on X ####
    dir_i = 0
    dirCase = 'X'
    # print(stripe_folder + ' ' + dirCase + ': start')
    stripe_folder_path = posixpath.join(msa_folders[dir_i], stripe_folder)
    # try:
    resultsX = collect_gmset_response(stripe_folder_path, beam_list_x, fracElement, dir_i + 1, spliceElement,
                                      splice_list_x, column_list_x, minrdrift=minrdrift,
                                      splice_frac_strain=splice_frac_strain)
    # except:
    #     print('ERROR: ' + model_name + '_' + stripe_folder + ' ' + dirCase)

    #### collect results on Y ####
    dir_i = 1
    dirCase = 'Y'
    # print(stripe_folder + ' ' + dirCase + ': start')
    stripe_folder_path = posixpath.join(msa_folders[dir_i], stripe_folder)
    # try:
    resultsY = collect_gmset_response(stripe_folder_path, beam_list_y, fracElement, dir_i + 1, spliceElement,
                                      splice_list_y, column_list_y, minrdrift=minrdrift,
                                      splice_frac_strain=splice_frac_strain)
    # except:
    #     print('ERROR: ' + model_name + '_' + stripe_folder + ' ' + dirCase)

    #### consolidate in a unique table ####
    # Get RSN for each ground motion
    record_namesX = resultsX.index
    record_namesY = resultsY.index
    rsnX = []
    for i in range(len(record_namesX)):
        rsnX.append(record_namesX[i].split('_')[0])
    rsnX = np.array(rsnX)
    rsnY = []
    for i in range(len(record_namesY)):
        rsnY.append(record_namesY[i].split('_')[0])
    rsnY = np.array(rsnY)

    # Replace dataframe indexes by RSN
    resultsX.index = rsnX
    resultsY.index = rsnY

    # Join the tables to include X and Y results
    results = pd.concat([resultsX, resultsY], axis=1, join="inner")

    # Merge EndCriteria columns
    endCriteria = results['EndCriteria'].values
    endCriteriaTotal = []
    for i in range(len(endCriteria)):
        if endCriteria[i, 0] == 'MaxDrift' or endCriteria[i, 1] == 'MaxDrift':
            endCriteriaTotal.append('MaxDrift')
        elif endCriteria[i, 0] == 'Inconvergence' or endCriteria[i, 1] == 'Inconvergence':
            endCriteriaTotal.append('Inconvergence')
        else:
            endCriteriaTotal.append('nonCollapse')
    results = results.drop('EndCriteria', axis=1)  # drop EndCriteria per direction
    results.insert(0, 'EndCriteriaX', endCriteria[:, 0])  # Insert EndCriteriaX with proper column name
    results.insert(0, 'EndCriteriaY', endCriteria[:, 1])  # Insert EndCriteriaY with proper column name
    results.insert(0, 'EndCriteria', endCriteriaTotal)  # Insert combined EndCriteria

    results.to_csv(results_filename)


def collect_single_response(model_name_all, stripe_folder_all, save_results_folder_all, msa_folder_all, beam_list_x_all,
                            beam_list_y_all, fracElement, spliceElement_all, splice_list_x_all, splice_list_y_all,
                            column_list_x_all,
                            column_list_y_all, minrdrift, splice_frac_strain, drift_out, rdrift_out, case_i):
    # INPUTS
    # Collects the EDP results considering each stripe of each case as an independent job
    #    model_name_all           = list of str with the case name to collect results from
    #    stripe_folder_all        = list of str with the foldername of the stripe to collect results from
    #    save_results_folder_all  = list of str with path to save results
    #    msa_folder_all          = list of list [path_to_results_X, path_to_results_Y]
    #    beam_list_x_all          = list of 2D np.array indicating which beams exist in the X frame
    #    beam_list_Y_all          = list of 2D np.array indicating which beams exist in the Y frame
    #    fracElement              = true  -> collect connections DS
    #                               false -> does NOT collect connections DS
    #
    #    spliceElement_all        = 1 -> collect splice data
    #                               0 -> do not collect splice data
    #    splice_list_x_all        = list of 2D np.array indicating the pier and story with splice in the X frame
    #    splice_list_y_all        = list of 2D np.array indicating the pier and story with splice in the Y frame
    #    column_list_x_all        = list of 2D np.array indicating which columns exist in the X frame
    #    column_list_y_all        = list of 2D np.array indicating which column exist in the X frame
    #    minrdrift                = minimum residual drift to collect
    #    splice_frac_strain       = strain limit for deciding fracture occured on splices
    #    drift_out          = Method to output peak drift
    #                         -> 'abs' = multiple columns with the peak absolute value for each floor
    #                         -> 'both' = multiple columns with the peak positive an negative value for each floor
    #    rdrift_out         = Method to output residual drift
    #                         -> 'max' = single column with the maximum residual across floors
    #                         -> 'all_abs' = multiple columns with the residual absolute value for each floor
    #                         -> 'all' = multiple columns with the residual for each floor with its sign
    #

    # Parse case to execute
    model_name = model_name_all[case_i]
    dirCase = model_name.split('dir')[1][0]
    stripe_folder = stripe_folder_all[case_i]
    save_results_folder = save_results_folder_all[case_i]
    msa_folder = msa_folder_all[case_i]

    if dirCase == 'X':
        spliceElement = spliceElement_all[case_i]
        beam_list = beam_list_x_all[case_i]
        column_list = column_list_x_all[case_i]
        splice_list = splice_list_x_all[case_i]
    else:
        spliceElement = spliceElement_all[case_i]
        beam_list = beam_list_y_all[case_i]
        column_list = column_list_y_all[case_i]
        splice_list = splice_list_y_all[case_i]

    # Define path to save results
    results_filename = posixpath.join(save_results_folder, 'EDP_' + model_name + '_' + stripe_folder + '.csv')

    #### collect results on singleDir ####
    stripe_folder_path = posixpath.join(msa_folder, stripe_folder)
    results = collect_gmset_response(stripe_folder_path, beam_list, fracElement, 1, spliceElement,
                                     splice_list, column_list, drift_out=drift_out, rdrift_out=rdrift_out, minrdrift=minrdrift,
                                     splice_frac_strain=splice_frac_strain)

    results.to_csv(results_filename)


def get_pz_response(results_folder, pz_list, filenames):
    # Read response for panel zones, currently takes the maximum of the time history
    #
    # INPUTS
    #    results_folder = path to folder with the results of NLRHA
    #    pz_list      = 2D np.array with 1 or 0 for the pz that exist
    #    filenames      = list with any group of the following alternatives
    #                     'pz_rot'
    #                     'all_disp': include time vector
    #
    # OUTPUTS
    #    pz_results = dictionary with all results for the panel zones
    #

    pz_results = dict(keys=filenames)
    n_stories, n_pier = pz_list.shape
    num_pz = np.sum(pz_list)

    # Read results as 1d array
    for file in filenames:
        filepath = posixpath.join(results_folder, file + '.out')
        # Read requested data (always at end-2 of analysis to take the last full result in case of inconvergence)
        with open(filepath) as f:
            last_line = 0
            prev_line = 0
            pp_line = 0
            for line in f:
                ppp_line = pp_line
                pp_line = prev_line
                prev_line = last_line
                if line == '\n' or type(line) == int: # skip blank lines
                    last_line = last_line
                else:
                    last_line = line
            aux = last_line.split()  # read second to last line in a time history since sometimes has errors
            try:
                aux = [float(x) for x in aux]  # convert to list of floats
            except:
                aux = np.array([0, 0])
            # only keeps set of results that include all connections (considering last 4 lines as options)
            if (file == 'pz_rot' and len(aux) < num_pz) or \
                    (file == 'all_disp' and len(aux) < 1 + num_pz):
                try:
                    aux = prev_line.split()
                    aux = [float(x) for x in aux]  # convert to list of floats
                except:
                    aux = np.array([0, 0])
                if (file == 'pz_rot' and len(aux) < num_pz) or \
                        (file == 'all_disp' and len(aux) < 1 + num_pz):
                    try:
                        aux = pp_line.split()
                        aux = [float(x) for x in aux]  # convert to list of floats
                    except:
                        aux = np.array([0, 0])
                    if (file == 'pz_rot' and len(aux) < num_pz) or \
                        (file == 'all_disp' and len(aux) < 1 + num_pz):
                        try:
                            aux = ppp_line.split()
                            aux = [float(x) for x in aux]  # convert to list of floats
                        except:
                            aux = np.array([0, 0])

            if file == 'all_disp':
                aux = aux[1:]
                aux = aux[0:num_pz]
                results_1d = np.abs(aux)
            if file == 'pz_rot':
                aux = aux[0:num_pz]
                results_1d = np.abs(aux)

        # If still not complete row of results fills the missing values with zeros
        # This work around is not bad because this occurs on collapsed cases
        if (file == 'pz_rot' and len(results_1d) < num_pz) or \
                    (file == 'all_disp' and len(results_1d) < num_pz):
            print('WARNING:Autocompleted PZ results: ' + results_folder)
            print('length of vector = ' + str(len(aux)) + '; should be = ' + str(num_pz))
            aux2 = np.zeros(num_pz)
            aux2[0:len(results_1d)] = results_1d
            results_1d = aux2

        # Initialize final matrix to return the results
        pz_results[file] = np.zeros([n_stories, n_pier])

        # Save in desired format
        i_element = 0
        for i_story in range(n_stories):
            for i_pier in range(n_pier):
                if pz_list[i_story, i_pier] > 0:
                    pz_results[file][i_story, i_pier] = results_1d[i_element]
                    i_element += 1

    return pz_results


def get_column_response(results_folder, column_list, filenames, def_desired='rot'):
    # Read response for columns, currently takes the maximum of the time history
    #
    # INPUTS
    #    results_folder = path to folder with the results of NLRHA
    #    beam_list      = 2D np.array with 1 or 0 for the beams that exist
    #    column_list    = 2D np.array with 1 or 0 for the columns that exist
    #    filenames      = list with any group of the following alternatives
    #                     'hinge_bot'
    #                     'hinge_top'
    #   def_desired     = 'axial'
    #                     'shear'
    #                     'rot'
    #                     Only applies for hinge_left/right and def_left/right
    #
    # OUTPUTS
    #    column_results = dictionary with all results for the columns, one key for each filename
    #

    column_results = dict(keys=filenames)

    n_stories, n_pier = column_list.shape
    num_columns = np.sum(column_list)

    # Read results as 1d array
    for file in filenames:
        filepath = posixpath.join(results_folder, file + '.out')

        # Read requested data (always at end-2 of analysis)
        with open(filepath) as f:
            last_line = 0
            prev_line = 0
            pp_line = 0
            for line in f:
                ppp_line = pp_line
                pp_line = prev_line
                prev_line = last_line
                if line == '\n'or type(line) == int:  # skip blank lines
                    last_line = last_line
                else:
                    last_line = line
            # try:
            aux = last_line.split()  # read last line in a time history
            # except:
                # print('ERROR: ' + filepath)

            try:
                aux = [float(x) for x in aux]  # convert to list of floats
            except:
                aux = np.array([0, 0])
            # only keeps lines with results for all the columns (possible issues when inconvergence)
            if ((len(aux) != int(3 * num_columns)) or (len(aux) != int(6 * num_columns)) or
                    (len(aux) != int(num_columns))):
                try:
                    aux = prev_line.split()
                    aux = [float(x) for x in aux]  # convert to list of floats
                except:
                    aux = np.array([0, 0])
                if ((len(aux) != int(3 * num_columns)) or (len(aux) != int(6 * num_columns)) or
                        (len(aux) != int(num_columns))):
                    try:
                        aux = pp_line.split()
                        aux = [float(x) for x in aux]  # convert to list of floats
                    except:
                        aux = np.array([0, 0])
                    if ((len(aux) != int(3 * num_columns)) or (len(aux) != int(6 * num_columns)) or
                            (len(aux) != int(num_columns))):
                        try:
                            aux = ppp_line.split()  # read second to last line in a time history since sometimes has errors
                            aux = [float(x) for x in aux]  # convert to list of floats
                        except:
                            aux = np.array([0, 0])
            data_1d = np.abs(aux)

        # If still not complete row of results fills the missing values with zeros
        # This work around is not bad because this occurs on collapsed cases
        if len(data_1d) < int(num_columns):
            print('WARNING:Autocompleted ' + file + ' results: ' + results_folder)
            print('length of vector = ' + str(len(data_1d)) + '; should be = ' + str(num_columns))
            aux2 = np.zeros(num_columns)
            aux2[0:len(data_1d)] = data_1d
            data_1d = aux2
        elif (len(data_1d) < int(3 * num_columns)) and (len(data_1d) != int(num_columns)):
            print('WARNING:Autocompleted ' + file + ' results: ' + results_folder)
            print('length of vector = ' + str(len(data_1d)) + '; should be = ' + str(3*num_columns))
            aux2 = np.zeros(3*num_columns)
            aux2[0:len(data_1d)] = data_1d
            data_1d = aux2
        elif (len(data_1d) < int(6 * num_columns)) and (len(data_1d) != int(num_columns)) and \
                (len(data_1d) != int(num_columns)):
            print('WARNING:Autocompleted ' + file + ' results: ' + results_folder)
            print('length of vector = ' + str(len(data_1d)) + '; should be = ' + str(6*num_columns))
            aux2 = np.zeros(6 * num_columns)
            aux2[0:len(data_1d)] = data_1d
            data_1d = aux2

        n_cols = len(data_1d)

        if n_cols == int(3 * num_columns):
            # read axial def, shear def, rotation for each hinge
            axial_def = data_1d[0:n_cols:3]
            shear_def = data_1d[1:n_cols:3]
            rot = data_1d[2:n_cols:3]

            if def_desired == 'axial':
                results_1d = axial_def
            elif def_desired == 'shear':
                results_1d = shear_def
            else:
                results_1d = rot

        elif n_cols == int(6 * num_columns):
            # read moment for each hinge
            M_bot = data_1d[2:n_cols:6]
            M_top = data_1d[5:n_cols:6]

            # read response at index t for each hinge
            if 'bot' in file:
                results_1d = abs(M_bot)
            elif 'top' in file:
                results_1d = abs(M_top)

        elif n_cols == num_columns:
            results_1d = data_1d

        else:
            print('ERROR: output not consistent with number of columns')
            print(results_folder)
            print(file)
            print('n_cols = ' + str(n_cols) + '; num_columns = ' + str(num_columns))

        # Initialize final matrix to return the results
        column_results[file] = np.zeros([n_stories, n_pier])

        # Save in desired format
        i_element = 0
        for i_story in range(n_stories):
            for i_pier in range(n_pier):
                if column_list[i_story, i_pier] > 0:  # jump if element does not exist in model setbacks
                    # if i_story == 0 or beam_list[
                    #     i_story - 1, min(i_pier, n_pier - 2)]:  # jump columns already created in atriums (no need since the column_list vector already includes this info)

                        column_results[file][i_story, i_pier] = results_1d[i_element]
                        i_element += 1

    return column_results


def get_beam_response(results_folder, beam_list, filenames, def_desired='rot'):
    # Read response for beams, currently takes either the maximum or the last of the time history
    #
    # INPUTS
    #    results_folder = path to folder with the results of NLRHA
    #    beam_list      = 2D np.array with 1 or 0 for the beams that exist
    #    filenames      = list with any group of the following alternatives
    #                     'frac_LB'
    #                     'frac_LT'
    #                     'frac_RB'
    #                     'frac_RT'
    #                     'FI_LB'
    #                     'FI_LT'
    #                     'FI_RB'
    #                     'FI_RT'
    #                     'hinge_left': plastic hinge on left end of beam
    #                     'hinge_right': plastic hinge on right end of beam
    #                     'def_left': fracture spring on left end of beam
    #                     'def_right': fracture spring on right end of beam
    #   def_desired     = 'axial'
    #                     'shear'
    #                     'rot'
    #                     Only applies for hinge_left/right and def_left/right
    #
    # OUTPUTS
    #    beam_results = dictionary with all results for the beams, one key for each filename
    #

    beam_results = dict(keys=filenames)

    n_stories, n_bays = beam_list.shape
    num_beams = np.sum(beam_list)

    # Read results as 1d array
    for file in filenames:
        filepath = posixpath.join(results_folder, file + '.out')

        # Read requested data (always at end-2 of analysis)
        with open(filepath) as f:
            last_line = 0
            prev_line = 0
            pp_line = 0
            for line in f:
                ppp_line = pp_line
                pp_line = prev_line
                prev_line = last_line
                if line == '\n' or type(line) == int:  # skip blank lines
                    last_line = last_line
                else:
                    last_line = line

            # try:
            aux = last_line.split()  # read last line in a time history
            # except:
            #     print('ERROR: ' + filepath)

            try:
                aux = [float(x) for x in aux]  # convert to list of floats
            except:
                aux = np.array([0, 0])
            # only keeps lines with results for all the columns (possible issues when inconvergence)
            if (len(aux) != int(3 * num_beams)) or (len(aux) != int(num_beams)):
                try:
                    aux = prev_line.split()
                    aux = [float(x) for x in aux]  # convert to list of floats
                except:
                    aux = np.array([0, 0])
                if (len(aux) != int(3 * num_beams)) or (len(aux) != int(num_beams)):
                    try:
                        aux = pp_line.split()
                        aux = [float(x) for x in aux]  # convert to list of floats
                    except:
                        aux = np.array([0, 0])
                    if (len(aux) != int(3 * num_beams)) or (len(aux) != int(num_beams)):
                        try:
                            aux = ppp_line.split()
                            aux = [float(x) for x in aux]  # convert to list of floats
                        except:
                            aux = np.array([0, 0])
            aux = [float(x) for x in aux]  # convert to list of floats
            data_1d = np.abs(aux)

        n_cols = len(data_1d)

        if n_cols == 3 * num_beams:  # file == 'hinge_left' or file == 'hinge_right':
            # read axial def, shear def, rotation for each hinge
            axial_def = data_1d[0:n_cols:3]
            shear_def = data_1d[1:n_cols:3]
            rot = data_1d[2:n_cols:3]

            if def_desired == 'axial':
                results_1d = axial_def
            elif def_desired == 'shear':
                results_1d = shear_def
            else:
                results_1d = rot

        elif n_cols == num_beams:
            results_1d = data_1d
        else:
            print('ERROR: output not consistent with number of beams')
            print(file)
            print('n_cols = ' + str(n_cols) + '; num_beams = ' + str(num_beams))
            beam_results = np.ones(num_beams) * np.nan
            return beam_results

        # Initialize final matrix to return the results
        beam_results[file] = np.zeros([n_stories, n_bays])

        # Save in desired format
        i_element = 0
        for i_story in range(n_stories):
            for i_beam in range(n_bays):
                if beam_list[i_story, i_beam] > 0:
                    beam_results[file][i_story, i_beam] = results_1d[i_element]
                    i_element += 1

    return beam_results


def get_splice_response(results_folder, splice_list, column_list, filenames, res_type='Max',
                        def_desired='rot'):
    # Read response for columns, currently takes the maximum of the time history
    #
    # INPUTS
    #    results_folder   = path to folder with the results of NLRHA
    #    splice_list      = 2D np.array with 1 or 0 for the beams that exist
    #    beam_list        = 2D np.array with 1 or 0 for the beams that exist
    #    column_list      = 2D np.array with 1 or 0 for the columns that exist
    #    filenames        = list with any group of the following alternatives
    #                       'ss_splice'
    #    res_type         = 'Max' : return the maximum response in the time history
    #                       'at_end': return the response at end of the history
    #    t                = index for deformed shape plot
    #   def_desired       = 'axial'
    #                       'shear'
    #                       'rot'
    #                       'strain'
    #                       'stress'
    #                       Only applies for hinge_left/right and def_left/right
    #
    # OUTPUTS
    #    column_results   = dictionary with all results for the columns, one key for each filename
    #

    column_results = dict(keys=filenames)

    n_stories, n_pier = column_list.shape
    # splice_list = np.zeros([n_stories, n_pier])
    # splice_list[n_stories_splice-1 : n_stories:n_stories_splice] = np.ones([n_pier])
    num_splices = np.sum(splice_list)
    max_res = 0

    # Read results as 1d array
    for file in filenames:
        filepath = posixpath.join(results_folder, file + '.out')

        # data_1d = np.loadtxt(filepath)
        # _, n_cols = data_1d.shape
        with open(filepath) as f:
            last_line = 0
            for line in f:
                prev_line = last_line
                try:
                    aux = line.split()
                    aux = [abs(float(x)) for x in aux]  # convert to list of floats, empty "line" gives an error here

                    # Initialize max_res vector if not exist yet
                    if type(max_res) == int:
                        max_res = np.zeros(len(aux))

                    # only keep row with complete data
                    if ((len(aux) == int(2 * num_splices)) or (len(aux) == int(10 * num_splices)) or
                            (len(aux) == int(3 * num_splices)) or (len(aux) == int(6 * num_splices)) or
                            (len(aux) == int(num_splices))):
                        last_line = line
                        max_res = np.max([max_res, aux], axis=0)
                    else:
                        last_line = last_line
                except:
                    last_line = last_line

            if res_type == 'Max':
                # Take the max value in the history
                data_1d = max_res
            else:
                # Taking last value
                aux = prev_line.split()
                aux = [float(x) for x in aux]  # convert to list of floats
                data_1d = np.abs(aux)

        n_cols = len(data_1d)

        if n_cols == int(2 * num_splices):
            # read stress and strain for a unique section
            stress = data_1d[0:n_cols:2]
            strain = data_1d[1:n_cols:2]

            if def_desired == 'stress':
                res = stress
            elif def_desired == 'strain':
                res = strain
            else:
                print('ERROR: specify either "stress" or "strain" output')

            results_1d = res

        elif n_cols == int(10 * num_splices):
            # read stress and strain for a unique section from results on all (5) section of the forceBasedElement
            stress = data_1d[:, 4:n_cols:10]
            strain = data_1d[:, 5:n_cols:10]

            if def_desired == 'stress':
                res = stress
            elif def_desired == 'strain':
                res = strain
            else:
                print('ERROR: specify either "stress" or "strain" output')

            results_1d = res

        elif n_cols == int(3 * num_splices):
            # read axial def, shear def, rotation for each hinge
            axial_def = data_1d[:, 0:n_cols:3]
            shear_def = data_1d[:, 1:n_cols:3]
            rot = data_1d[:, 2:n_cols:3]

            if def_desired == 'axial':
                res = axial_def
            elif def_desired == 'shear':
                res = shear_def
            else:
                res = rot

            results_1d = res

        elif n_cols == int(6 * num_splices):
            # read moment for each hinge
            M_bot = data_1d[:, 2:n_cols:6]
            M_top = data_1d[:, 5:n_cols:6]

            # read response history
            if 'bot' in file:
                aux = M_bot
            elif 'top' in file:
                aux = M_top

        elif n_cols == num_splices:
            results_1d = data_1d
        else:
            print('ERROR: output not consistent with number of columns')
            print(file)
            print('n_cols = ' + str(n_cols) + '; num_splices = ' + str(num_splices))

        # Initialize final matrix to return the results
        if res_type == 'all_t':
            n_pts, _ = results_1d.shape
            column_results[file] = np.zeros([n_stories, n_pier, n_pts])
        else:
            column_results[file] = np.zeros([n_stories, n_pier])

        # Save in desired format
        # i_element = 0
        # for i_pier in range(n_pier):
        #     for i_story in range(n_stories):
        #         if column_list[i_story, i_pier] > 0:  # jump if setbacks
        #             if i_story == 0 or beam_list[
        #                 i_story - 1, min(i_pier, n_pier - 2)]:  # jump columns already created in atriums
        #                 if splice_list[i_story, i_pier] > 0:  # jump if no splice in this story
        #                     if res_type == 'all_t':
        #                         column_results[file][i_story, i_pier, :] = results_1d[:, i_element]
        #                     else:
        #                         column_results[file][i_story, i_pier] = results_1d[i_element]
        #
        #                     i_element += 1
        i_element = 0
        for i_pier in range(n_pier):
            for i_story in range(n_stories):
                if splice_list[i_story, i_pier] > 0:
                # jump if no splice in this column segment. The splice_list array already accounts for setbacks, atriums or missing columns
                    if res_type == 'all_t':
                        column_results[file][i_story, i_pier, :] = results_1d[:, i_element]
                    else:
                        column_results[file][i_story, i_pier] = results_1d[i_element]

                    i_element += 1

    return column_results


def collect_endState_single_response(model_name_all, save_results_folder_all, stripe_folders_all, msa_folder_all,
                                     beam_list_x_all,
                                     beam_list_y_all, column_list_x_all, column_list_y_all, pz_list_x_all,
                                     pz_list_y_all, splice_all, colSplice_x_all,
                                     colSplice_y_all, case_i):
    # INPUTS
    # All the inputs include information per case (different from the EDP collector that breaks each case into independent jobs per stripe
    #    model_name_all           = list of str with the case name to collect results from
    #    save_results_folder_all  = list of str with path to save results
    #    stripe_folders_all       = list of list with the folder name of each stripe
    #    msa_folder_all          = list of list [path_to_results_X, path_to_results_Y]
    #    beam_list_x_all          = list of 2D np.array indicating which beams exist in the X frame
    #    beam_list_Y_all          = list of 2D np.array indicating which beams exist in the Y frame
    #    column_list_x_all        = list of 2D np.array indicating which columns exist in the X frame
    #    column_list_y_all        = list of 2D np.array indicating which columns exist in the Y frame
    #    pz_list_x_all            = list of 2D np.array indicating which pz exist in the X frame
    #    pz_list_y_all            = list of 2D np.array indicating which pz exist in the Y frame
    #    splice_all               = list of boolean if splice are considered or ignored
    #    colSplice_x_all          = list of 2D np.array indicating which stories have a splice in the X frame
    #    colSplice_y_all          = list of 2D np.array indicating which stories have a splice in the X frame
    #

    # Parse case to execute
    model_name = model_name_all[case_i]
    dirCase = model_name.split('dir')[1][0]
    save_results_folder = save_results_folder_all[case_i]
    stripe_folders = stripe_folders_all[case_i]
    msa_folder = msa_folder_all[case_i]
    beam_list_x = beam_list_x_all[case_i]
    beam_list_y = beam_list_y_all[case_i]
    column_list_x = column_list_x_all[case_i]
    column_list_y = column_list_y_all[case_i]
    pz_list_x = pz_list_x_all[case_i]
    pz_list_y = pz_list_y_all[case_i]
    splice = splice_all[case_i]
    if splice == 1:
        colSplice_x = colSplice_x_all[case_i]
        colSplice_y = colSplice_y_all[case_i]

    n_stripes = len(stripe_folders)

    # Collect end state for every gm in every return period for GIVEN DIRECTIONS
    print('------- ' + model_name + ' FRAME IN ' + dirCase + ' -------')

    # Select frame geometry matrices
    if dirCase == 'X':
        beam_list = beam_list_x
        column_list = column_list_x
        pz_list = pz_list_x
        if splice == 1:
            colSplice = colSplice_x
    else:
        beam_list = beam_list_y
        column_list = column_list_y
        pz_list = pz_list_y
        if splice == 1:
            colSplice = colSplice_y
    results_filename = posixpath.join(save_results_folder, 'end_state_' + model_name + '_dir' + dirCase + '.h5')

    # # count panel zones
    n_stories, n_pier = column_list.shape
    # num_pz = 0
    # for i_story in range(n_stories):
    #     for i_pier in range(n_pier):
    #         if ((column_list[min(i_story + 1, n_stories-1), i_pier] == 1) or
    #                 ((column_list[i_story, i_pier] == 1 )) and
    #                 (beam_list[i_story, i_pier - 1])):
    #             existPZ = True
    #         elif ((column_list[min(i_story + 1, n_stories-1), i_pier] == 1) or
    #                 ((column_list[i_story, i_pier] == 1 ) and
    #                  (beam_list[i_story, min(i_pier, n_pier -2)]))):
    #             existPZ = True
    #         else:
    #             existPZ = False
    #
    #         if existPZ:
    #             num_pz += 1
    # print('num_pz='+str(num_pz))

    # Removes existing file
    # if os.path.isfile(results_filename):
    #     os.remove(results_filename)
    #     print(results_filename + ' already exists, so deleted it')

    # if True: (Collects only data for those building without data)
    if not os.path.isfile(results_filename):

        # Collect results and store in HDF file
        with h5py.File(results_filename, 'w') as hf:
            # prepare data groups per return period
            for group in stripe_folders:
                _ = hf.create_group('/' + group)

            # collect bldg response for each gm in each stripe
            for i in range(n_stripes):
                stripe_folder_path = posixpath.join(msa_folder, stripe_folders[i])

                print('RP = ' + str(stripe_folders[i]) + 'years')
                # print(stripe_folder_path)

                gm_ids = os.listdir(stripe_folder_path)
                for j in range(len(gm_ids)):
                    # print(gm_ids[j])

                    # collect results for this gm
                    results_gm = posixpath.join(stripe_folder_path, gm_ids[j])

                    # check if acc results available (gm finished?)
                    pfa_gm = get_EDPstory_response(results_gm, n_stories, 'acc_env')
                    if type(pfa_gm) == int:  # did not finish RHA, so skip the ground motion
                        print('Did not finish GM' + str(gm_ids[j]))
                    else:
                        #    Panel zones
                        # pz_response     = get_pz_response(results_gm, beam_list, column_list, num_pz, ['all_disp', 'pz_rot'])
                        pz_response = get_pz_response(results_gm, pz_list, ['all_disp', 'pz_rot'])
                        #    beams and columns
                        if 'cvn' in model_name:
                            column_response = get_column_response(results_gm, column_list, ['hinge_bot', 'hinge_top'])
                            beam_plas_rot = get_beam_response(results_gm, beam_list, ['hinge_left', 'hinge_right'])
                            frac_simulated = get_beam_response(results_gm, beam_list,
                                                               ['frac_LB', 'frac_LT', 'frac_RB', 'frac_RT'])
                            FI_simulated = get_beam_response(results_gm, beam_list,
                                                               ['FI_LB', 'FI_LT', 'FI_RB', 'FI_RT'])
                        else:
                            beam_plas_rot = get_beam_response(results_gm, beam_list, ['hinge_left', 'hinge_right'])
                            column_response = get_column_response(results_gm, column_list, ['hinge_bot', 'hinge_top'])
                        #    Splices
                        if splice == 1:
                            splice_response = get_splice_response(results_gm, colSplice, column_list, ['ss_splice'],
                                                                  res_type='Max', def_desired='strain')
                            splice_frac = splice_response['ss_splice'] > 2 * 60 / 29000

                        # create gm group
                        rp_group = hf['/' + stripe_folders[i]]
                        gm_record_group = rp_group.create_group(gm_ids[j])

                        # Save in h5 file's building_group
                        if 'cvn' in model_name:
                            key = 'all_disp'
                            _ = gm_record_group.create_dataset(key, data=pz_response[key])
                            key = 'pz_rot'
                            _ = gm_record_group.create_dataset(key, data=pz_response[key])
                            key = 'hinge_bot'
                            _ = gm_record_group.create_dataset(key, data=column_response[key])
                            key = 'hinge_top'
                            _ = gm_record_group.create_dataset(key, data=column_response[key])
                            key = 'hinge_left'
                            _ = gm_record_group.create_dataset(key, data=beam_plas_rot[key])
                            key = 'hinge_right'
                            _ = gm_record_group.create_dataset(key, data=beam_plas_rot[key])
                            key = 'frac_LB'
                            _ = gm_record_group.create_dataset(key, data=frac_simulated[key])
                            key = 'frac_LT'
                            _ = gm_record_group.create_dataset(key, data=frac_simulated[key])
                            key = 'frac_RB'
                            _ = gm_record_group.create_dataset(key, data=frac_simulated[key])
                            key = 'frac_RT'
                            _ = gm_record_group.create_dataset(key, data=frac_simulated[key])
                            key = 'FI_LB'
                            _ = gm_record_group.create_dataset(key, data=FI_simulated[key])
                            key = 'FI_LT'
                            _ = gm_record_group.create_dataset(key, data=FI_simulated[key])
                            key = 'FI_RB'
                            _ = gm_record_group.create_dataset(key, data=FI_simulated[key])
                            key = 'FI_RT'
                            _ = gm_record_group.create_dataset(key, data=FI_simulated[key])
                            if splice == 1:
                                key = 'ss_splice'
                                _ = gm_record_group.create_dataset(key, data=splice_response[key])
                                key = 'splice_frac'
                                _ = gm_record_group.create_dataset(key, data=splice_frac)
                        else:
                            key = 'all_disp'
                            _ = gm_record_group.create_dataset(key, data=pz_response[key])
                            key = 'pz_rot'
                            _ = gm_record_group.create_dataset(key, data=pz_response[key])
                            key = 'hinge_bot'
                            _ = gm_record_group.create_dataset(key, data=column_response[key])
                            key = 'hinge_top'
                            _ = gm_record_group.create_dataset(key, data=column_response[key])
                            key = 'hinge_left'
                            _ = gm_record_group.create_dataset(key, data=beam_plas_rot[key])
                            key = 'hinge_right'
                            _ = gm_record_group.create_dataset(key, data=beam_plas_rot[key])
                            if splice == 1:
                                key = 'ss_splice'
                                _ = gm_record_group.create_dataset(key, data=splice_response[key])
                                key = 'splice_frac'
                                _ = gm_record_group.create_dataset(key, data=splice_frac)


def collect_endStateXandY_response(model_name_all, save_results_folder_all, stripe_folders_all, msa_folders_all, beam_list_x_all,
                                   beam_list_y_all, column_list_x_all, column_list_y_all, pz_list_x_all, pz_list_y_all, splice_all, colSplice_x_all,
                                   colSplice_y_all, case_i):
    # INPUTS
    # All the inputs include information per case (different from the EDP collector that breaks each case into independent jobs per stripe
    #    model_name_all           = list of str with the case name to collect results from
    #    save_results_folder_all  = list of str with path to save results
    #    stripe_folders_all       = list of list with the folder name of each stripe
    #    msa_folders_all          = list of list [path_to_results_X, path_to_results_Y]
    #    beam_list_x_all          = list of 2D np.array indicating which beams exist in the X frame
    #    beam_list_Y_all          = list of 2D np.array indicating which beams exist in the Y frame
    #    column_list_x_all        = list of 2D np.array indicating which columns exist in the X frame
    #    column_list_y_all        = list of 2D np.array indicating which columns exist in the Y frame
    #    pz_list_x_all            = list of 2D np.array indicating which pz exist in the X frame
    #    pz_list_y_all            = list of 2D np.array indicating which pz exist in the Y frame
    #    splice_all               = list of boolean if splice are considered or ignored
    #    colSplice_x_all          = list of 2D np.array indicating which stories have a splice in the X frame
    #    colSplice_y_all          = list of 2D np.array indicating which stories have a splice in the X frame
    #

    # Parse case to execute
    model_name    = model_name_all[case_i]
    save_results_folder = save_results_folder_all[case_i]
    stripe_folders= stripe_folders_all[case_i]
    msa_folders   = msa_folders_all[case_i]
    beam_list_x   = beam_list_x_all[case_i]
    beam_list_y   = beam_list_y_all[case_i]
    column_list_x = column_list_x_all[case_i]
    column_list_y = column_list_y_all[case_i]
    pz_list_x     = pz_list_x_all[case_i]
    pz_list_y     = pz_list_y_all[case_i]
    splice        = splice_all[case_i]
    if splice == 1:
        colSplice_x = colSplice_x_all[case_i]
        colSplice_y = colSplice_y_all[case_i]

    dirCases = ['X', 'Y']
    n_stripes = len(stripe_folders)
    # Collect end state for every gm in every return period for BOTH DIRECTIONS
    for dir_i in range(len(dirCases)):

        print('------- ' + model_name + ' FRAME IN ' + dirCases[dir_i] + ' -------')

        # Select frame geometry matrices
        if dir_i == 0:
            beam_list   = beam_list_x
            column_list = column_list_x
            pz_list = pz_list_x
            if splice == 1:
                colSplice = colSplice_x
        else:
            beam_list   = beam_list_y
            column_list = column_list_y
            pz_list = pz_list_y
            if splice == 1:
                colSplice = colSplice_y
        results_filename = posixpath.join(save_results_folder, 'end_state_' + model_name + '_dir'+ dirCases[dir_i] +'.h5')

        # # count panel zones
        n_stories, n_pier = column_list.shape
        # num_pz = 0
        # for i_story in range(n_stories):
        #     for i_pier in range(n_pier):
        #         if ((column_list[min(i_story + 1, n_stories-1), i_pier] == 1) or
        #                 ((column_list[i_story, i_pier] == 1 )) and
        #                 (beam_list[i_story, i_pier - 1])):
        #             existPZ = True
        #         elif ((column_list[min(i_story + 1, n_stories-1), i_pier] == 1) or
        #                 ((column_list[i_story, i_pier] == 1 ) and
        #                  (beam_list[i_story, min(i_pier, n_pier -2)]))):
        #             existPZ = True
        #         else:
        #             existPZ = False
        #
        #         if existPZ:
        #             num_pz += 1
        # print('num_pz='+str(num_pz))

        # # Removes existing file
        # if os.path.isfile(results_filename):
        #     os.remove(results_filename)
        #     print(results_filename + ' already exists, so deleted it')

        # if True (Collects only data for those building without data)
        if not os.path.isfile(results_filename):

            # Collect results and store in HDF file
            with h5py.File(results_filename, 'w') as hf:
                # prepare data groups per return period
                for group in stripe_folders:
                    _ = hf.create_group('/' + group)

                # collect bldg response for each gm in each stripe
                for i in range(n_stripes):
                    stripe_folder_path = posixpath.join(msa_folders[dir_i], stripe_folders[i])

                    print('RP = ' + str(stripe_folders[i]) + 'years')
                    # print(stripe_folder_path)

                    gm_ids = os.listdir(stripe_folder_path)
                    for j in range(len(gm_ids)):
                        # print(gm_ids[j])

                        # collect results for this gm
                        results_gm = posixpath.join(stripe_folder_path, gm_ids[j])

                        # check if acc results available (gm finished?)
                        pfa_gm = get_EDPstory_response(results_gm, n_stories, 'acc_env')
                        if type(pfa_gm) == int:  # did not finish RHA, so skip the ground motion
                            print('Did not finish GM' + str(gm_ids[j]))
                        else:
                            #    Panel zones
                            # pz_response     = get_pz_response(results_gm, beam_list, column_list, num_pz, ['all_disp', 'pz_rot'])
                            pz_response = get_pz_response(results_gm, pz_list, ['all_disp', 'pz_rot'])
                            #    beams and columns
                            if 'cvn' in model_name:
                                column_response = get_column_response(results_gm, column_list, ['hinge_bot','hinge_top'])
                                beam_plas_rot = get_beam_response(results_gm, beam_list, ['hinge_left', 'hinge_right'])
                                frac_simulated  = get_beam_response(results_gm, beam_list, ['frac_LB','frac_LT','frac_RB','frac_RT'])
                            else:
                                beam_plas_rot   = get_beam_response(results_gm, beam_list, ['hinge_left','hinge_right'])
                                column_response = get_column_response(results_gm, column_list, ['hinge_bot','hinge_top'])
                            #    Splices
                            if splice == 1:
                                splice_response = get_splice_response(results_gm, colSplice, column_list, ['ss_splice'],
                                                  res_type='Max', def_desired='strain')
                                splice_frac = splice_response['ss_splice'] > 2*60/29000

                            # create gm group
                            rp_group = hf['/' + stripe_folders[i]]
                            gm_record_group = rp_group.create_group(gm_ids[j])

                            # Save in h5 file's building_group
                            if 'cvn' in model_name:
                                key = 'all_disp'
                                _ = gm_record_group.create_dataset(key, data=pz_response[key])
                                key = 'pz_rot'
                                _ = gm_record_group.create_dataset(key, data=pz_response[key])
                                key = 'hinge_bot'
                                _ = gm_record_group.create_dataset(key, data=column_response[key])
                                key = 'hinge_top'
                                _ = gm_record_group.create_dataset(key, data=column_response[key])
                                key = 'hinge_left'
                                _ = gm_record_group.create_dataset(key, data=beam_plas_rot[key])
                                key = 'hinge_right'
                                _ = gm_record_group.create_dataset(key, data=beam_plas_rot[key])
                                key = 'frac_LB'
                                _ = gm_record_group.create_dataset(key, data=frac_simulated[key])
                                key = 'frac_LT'
                                _ = gm_record_group.create_dataset(key, data=frac_simulated[key])
                                key = 'frac_RB'
                                _ = gm_record_group.create_dataset(key, data=frac_simulated[key])
                                key = 'frac_RT'
                                _ = gm_record_group.create_dataset(key, data=frac_simulated[key])
                                if splice == 1:
                                    key = 'ss_splice'
                                    _ = gm_record_group.create_dataset(key, data=splice_response[key])
                                    key = 'splice_frac'
                                    _ = gm_record_group.create_dataset(key, data=splice_frac)
                            else:
                                key = 'all_disp'
                                _ = gm_record_group.create_dataset(key, data=pz_response[key])
                                key = 'pz_rot'
                                _ = gm_record_group.create_dataset(key, data=pz_response[key])
                                key = 'hinge_bot'
                                _ = gm_record_group.create_dataset(key, data=column_response[key])
                                key = 'hinge_top'
                                _ = gm_record_group.create_dataset(key, data=column_response[key])
                                key = 'hinge_left'
                                _ = gm_record_group.create_dataset(key, data=beam_plas_rot[key])
                                key = 'hinge_right'
                                _ = gm_record_group.create_dataset(key, data=beam_plas_rot[key])
                                if splice == 1:
                                    key = 'ss_splice'
                                    _ = gm_record_group.create_dataset(key, data=splice_response[key])
                                    key = 'splice_frac'
                                    _ = gm_record_group.create_dataset(key, data=splice_frac)
