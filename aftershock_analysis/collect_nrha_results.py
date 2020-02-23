from .base import *


def store_building_geometry(results_filename, hf_group, n_stories, n_bays, story_height, bay_width):

    with h5py.File(results_filename, 'r+') as hf:
        hf_group = hf[hf_group]

        hf_group.attrs['n_stories'] = n_stories
        hf_group.attrs['n_floors'] = n_stories + 1

        hf_group.attrs['n_bays'] = n_bays
        hf_group.attrs['n_columns'] = n_bays + 1

        bay_width = bay_width * 12
        story_height = story_height * 12

        hf_group.attrs['bay_width'] = bay_width
        hf_group.attrs['story_height'] = n_stories * story_height

        hf_group.attrs['building_width'] = n_bays * bay_width
        hf_group.attrs['building_height'] = n_stories * story_height

        # store the original geometry of each column
        columns = np.zeros(((n_bays + 1) * n_stories, 2, 2))
        i_element = 0
        for i_story in range(n_stories):
            for i_beam in range(n_bays + 1):
                # x values of columns
                columns[i_element, :, 0] = (i_beam) * bay_width
                for i_end in range(2):
                    # y values of columns
                    columns[i_element, i_end, 1] = (i_story + i_end) * story_height
                i_element = i_element + 1
        key = 'column_geometry'
        hf_group.create_dataset(key, data=columns)

        # store the original geometry of each beam
        beams = np.zeros((n_bays * n_stories, 2, 2))
        i_element = 0
        for i_story in range(n_stories):
            for i_beam in range(n_bays):
                # y values of beams
                beams[i_element, :, 1] = (i_story + 1) * story_height
                for i_end in range(2):
                    # x values of beams
                    beams[i_element, i_end, 0] = (i_beam + i_end) * bay_width
                i_element = i_element + 1
        key = 'beam_geometry'
        hf_group.create_dataset(key, data=beams)


def collect_gm_metadata(gm_files, results_filename, gm_group):

    [gm_metadata_file,
     gm_sa_avg_file,
     gm_sa_t1_file,
     gm_duration_file,
     gm_spectra_file,
     gm_acc_folder] = gm_files

    gm_metadata = pd.read_csv(gm_metadata_file, sep='\t')
    n_gms = len(gm_metadata)
    gm_ids = ['GM' + str(i + 1) for i in range(n_gms)]
    gm_metadata['id'] = gm_ids
    gm_metadata.set_index('id', inplace=True)

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
    gm_metadata['Duration_5-75'] = durations[' t_s575 '].values

    response_spectra = pd.read_csv(gm_spectra_file, sep='\t', index_col='Record ')
    new_col = [x.strip() for x in response_spectra.columns]
    new_index = [x.strip() for x in response_spectra.index]
    response_spectra['id'] = new_index
    response_spectra.set_index('id', inplace=True)
    response_spectra.columns = new_col
    periods = [float(x[3:-2]) for x in new_col]

    key = '/ground_motion_records/gm_response_spectra'
    response_spectra.to_hdf(results_filename, key=key)


    with h5py.File(results_filename, 'r+') as hf:
        gm_group = hf[gm_group]
        gm_group['gm_response_spectra'].attrs['periods'] = periods

        # store each ground motion
        for gm_id in gm_ids:
            gm_record_group = gm_group.create_group(gm_id)

            gm_record_group.attrs['rsn'] = gm_metadata.loc[gm_id, 'RSN']
            gm_record_group.attrs['event'] = gm_metadata.loc[gm_id, 'eventName']
            gm_record_group.attrs['date'] = gm_metadata.loc[gm_id, 'Date']
            gm_record_group.attrs['station'] = gm_metadata.loc[gm_id, 'Station']
            gm_record_group.attrs['magnitude'] = gm_metadata.loc[gm_id, 'M']
            gm_record_group.attrs['r_rup'] = gm_metadata.loc[gm_id, 'Rup']
            gm_record_group.attrs['r_jb'] = gm_metadata.loc[gm_id, 'Rjb']
            gm_record_group.attrs['vs30'] = gm_metadata.loc[gm_id, 'Vs30']
            gm_record_group.attrs['region'] = gm_metadata.loc[gm_id, 'region']
            gm_record_group.attrs['fault_type'] = gm_metadata.loc[gm_id, 'Fault_Type']
            gm_record_group.attrs['component'] = gm_metadata.loc[gm_id, 'Component']
            gm_record_group.attrs['unscaled_sa_t1'] = gm_metadata.loc[gm_id, 'Unscaled Sa(T1)']
            gm_record_group.attrs['unscaled_sa_avg'] = gm_metadata.loc[gm_id, 'Unscaled Sa_avg']

            # acceleration time history
            gm_acc_file = posixpath.join(gm_acc_folder, gm_id + '.txt')
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


def collect_ida_results(ida_folder, gm_metadata, results_filename, ida_results_group):
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


def collect_damaged_results(damaged_folder, gm_metadata, results_filename, damaged_group, building_group, result_type):
    gm_ids = gm_metadata.index

    if result_type != 'mainshock_edp':

        all_mainshocks = os.listdir(damaged_folder)

        for gm_id in gm_ids:

            gm_id_mainshocks = [i for i in all_mainshocks if gm_id + '_' in i]
            scales = np.sort([float(x[-6:-3]) for x in gm_id_mainshocks])

            for scale in scales:
                gm_scale_group = create_damaged_gm_scale_group(results_filename, damaged_group, gm_id, scale, gm_metadata)

                scale_name = str(scale) + 'Col'
                gm_scale_folder = posixpath.join(damaged_folder, gm_id + '_' + scale_name)
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
                    collect_ida_time_history(gm_scale_folder, gm_metadata, results_filename, ida_results_group)
                else:
                    raise ValueError('Add code for result_type.')

    else:

        scales = np.sort([float(x[3:]) for x in os.listdir(damaged_folder)])

        remove_idx = scales==1.0
        scales = scales[~remove_idx]

        peak_drift_max = pd.DataFrame(index=gm_ids, columns=scales, dtype='float64')
        residual_drift_max = pd.DataFrame(index=gm_ids, columns=scales, dtype='float64')

        for scale in scales:
            scale_name = 'STR' + str(scale) + '0'

            for gm_id in gm_ids:
                gm_scale_group = create_damaged_gm_scale_group(results_filename, damaged_group, gm_id, scale, gm_metadata)
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

        scale_name = str(scale) + 'Col'
        if scale_name in gm_damaged_group.keys():
            gm_scale_group = gm_damaged_group[scale_name]
        else:
            gm_scale_group = gm_damaged_group.create_group(scale_name)

            collapse_intensity = gm_metadata.loc[gm_id, ['Intact Collapse Scale Factor',
                                                         'Intact Collapse Sa(T1)',
                                                         'Intact Collapse Sa_avg']].values

            gm_scale_group.attrs['scale_factor'] = scale * collapse_intensity[0]
            gm_scale_group.attrs['sa_t1'] = scale * collapse_intensity[1]
            gm_scale_group.attrs['sa_avg'] = scale * collapse_intensity[2]

        return gm_scale_group.name


def collect_mainshock_edp_results(edp_folder, hf, building_group, edp_results_group):
    edp_list = ['drift', 'displacement', 'acceleration']
    edp_list = ['drift', 'displacement']

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
            dset = edp_results_group.create_dataset(edp_name, data=edp_results)

            edp_name = 'peak_story_drift'
            peak_results = np.max(np.abs(edp_results), axis=1)
            dset = edp_results_group.create_dataset(edp_name, data=peak_results)
            dset.attrs['max_peak_story_drift'] = np.max(peak_results)

            edp_name = 'residual_story_drift'
            residual_results = np.median(edp_results[:, residual_start_idx:-1], axis=1)
            dset = edp_results_group.create_dataset(edp_name, data=residual_results)
            dset.attrs['max_residual_story_drift'] = np.max(np.abs(residual_results))

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
            dset = edp_results_group.create_dataset(edp_name, data=edp_results)
            dset.attrs['units'] = 'inches'

            edp_name = 'peak_displacement'
            peak_results = np.max(np.abs(edp_results), axis=1)
            dset = edp_results_group.create_dataset(edp_name, data=peak_results)
            dset.attrs['max_peak_displacement'] = np.max(peak_results)
            dset.attrs['units'] = 'inches'

            edp_name = 'residual_displacement'
            residual_results = np.median(edp_results[:, residual_start_idx:-1], axis=1)
            dset = edp_results_group.create_dataset(edp_name, data=residual_results)
            dset.attrs['max_residual_displacement'] = np.max(np.abs(residual_results))
            dset.attrs['units'] = 'inches'

        else:
            raise ValueError('define edp results collection method')