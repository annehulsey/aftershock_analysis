from .base import *
from .collect_nrha_results import *


def collect_ida_curves(results_filename, gm_metadata, ida_folder):
    gm_ids = gm_metadata.index
    n_gms = len(gm_ids)

    ida_segments = np.zeros((n_gms, 100, 2))
    ida_segments[:, :, 0] = 0.1

    for i, gm_id in zip(range(n_gms), gm_ids):
        key = ida_folder + gm_id + '/ida_curve'
        ida_curve = pd.read_hdf(results_filename, key)

        ida_curve = ida_curve.loc[:, ['Sa_avg', 'Story Drift Ratio (max)']].to_numpy()
        collapse = ida_curve[-1, 0]
        ida_segments[i, :, 1] = collapse

        ida_segments[i, 0:len(ida_curve), 1] = ida_curve[:, 0]
        ida_segments[i, 0:len(ida_curve), 0] = ida_curve[:, 1]

    return ida_segments


def collect_peak_and_residual_drift_curves(results_filename, gm_metadata):
    gm_ids = gm_metadata.index
    n_gms = len(gm_ids)
    collapse_intensities = gm_metadata['Intact Collapse Sa_avg']

    key = '/mainshock_damage_results/peak_story_drift_max'
    peak_drift = pd.read_hdf(results_filename, key=key)
    peak_drift = peak_drift[np.sort(peak_drift.columns)].copy()

    key = '/mainshock_damage_results/residual_drift_max'
    residual_drift = pd.read_hdf(results_filename, key=key)
    residual_drift = residual_drift[np.sort(peak_drift.columns)].copy()

    stripes = peak_drift.columns
    n_stripes = len(stripes)

    for edp_type in ['peak', 'residual']:
        if edp_type == 'peak':
            edp = peak_drift
        elif edp_type == 'residual':
            edp = residual_drift

        segments = np.zeros((n_gms, n_stripes + 1, 2))
        segments[:, 0, :] = 0
        segments[:, 1:, 0] = edp.to_numpy()

        for stripe, i in zip(stripes, range(n_stripes)):
            segments[:, i + 1, 1] = collapse_intensities * stripe

        if edp_type == 'peak':
            peak_segments = segments
        elif edp_type == 'residual':
            residual_segments = segments

    return peak_segments, residual_segments


def plot_damaged_msa_vs_ida_per_gm(gm_id, results_filename, gm_metadata, msa_intact_fragility, ida_intact_fragility, intact_ida_segments):

    damaged_group = 'mainshock_damage_results'
    fragility_linewidth = 3

    with h5py.File(results_filename, 'r') as hf:
        i = int(gm_id[2:]) - 1

        msa_fragilities = msa_intact_fragility.copy()
        ida_fragilities = ida_intact_fragility.copy()
        stripes = list(hf[damaged_group][gm_id].keys())

        if len(stripes) > 1:
            fig, ax = plt.subplots(2,2,figsize=(15,10), gridspec_kw={'width_ratios':[0.55,0.45]})
            suptitle = fig.suptitle('Mainshock: ' + gm_id)
            ax[0][1].get_shared_x_axes().join(ax[0][1], ax[1][1])

            # intact ida
            current_ax = ax[0][0]
            ida_ylim = [0, 1.1*np.amax(intact_ida_segments[:,:,1])]
            full_ida_plot = LineCollection(intact_ida_segments, linewidths=1, colors='lightgray', linestyle='solid')
            _ = current_ax.add_collection(full_ida_plot)
            _ = current_ax.set_xlim(0,0.1)
            _ = current_ax.set_ylim(ida_ylim)
            _ = current_ax.set_ylabel('Sa$_{avg}$(T$_1$), [g]')
            _ = current_ax.set_xlabel('Peak Story Drift Ratio')

            _ = current_ax.plot(intact_ida_segments[i,:,0], intact_ida_segments[i,:,1], color='k', linewidth=2.5, label=gm_id)
            _ = current_ax.legend(loc='upper left')

            if True:
                current_ax = ax[1][0]
                full_ida_plot = LineCollection(intact_ida_segments, linewidths=1, colors='lightgray', linestyle='solid')
                _ = current_ax.add_collection(full_ida_plot)


            # intact msa fragilities
            current_ax = ax[0][1]
            median = msa_fragilities.loc['Intact','Median']
            beta = msa_fragilities.loc['Intact','Beta']

            s = 0
            y = np.linspace(0.001, 1, num=100)
            x = stats.lognorm(beta, scale=median).ppf(y)
            label = 'No Mainshock'+ '\n$IM_{0.5}=$'+'{0:.2f}'.format(median) + ', $\sigma_{ln}=$'+'{0:.2f}'.format(beta)
            _ = current_ax.plot(x, 100*y, label=label, color='C'+str(s), linewidth=fragility_linewidth)

            key = '/intact_results/msa_sa_avg/collapse_matrix'
            msa_stripes = pd.read_hdf(results_filename, key).columns
            collapse_matrix = pd.read_hdf(results_filename, key).to_numpy()
            [n_gms, _] = collapse_matrix.shape
            percent_collapse = 100*np.sum(collapse_matrix, axis=0)/n_gms
            _ = current_ax.scatter(msa_stripes, percent_collapse, color='C'+str(s), edgecolors='k', linewidth=1)

            # intact ida fragilities
            current_ax = ax[1][1]
            median = ida_fragilities.loc['Intact','Median']
            beta = ida_fragilities.loc['Intact','Beta']

            y = np.linspace(0.001, 1, num=100)
            x = stats.lognorm(beta, scale=median).ppf(y)
            label = 'No Mainshock'+ '\n$IM_{0.5}=$'+'{0:.2f}'.format(median) + ' $\sigma_{ln}=$'+'{0:.2f}'.format(beta)
            _ = current_ax.plot(x, 100*y, label=label, color='C'+str(s), linewidth=2.5)

            key = '/intact_results/ida/collapse_intensities'
            collapse_intensities = pd.read_hdf(results_filename, key)['Sa_avg'].to_numpy()
            _ = current_ax.scatter(np.sort(collapse_intensities), np.linspace(100/n_gms, 100, num=n_gms, endpoint=True),
                                   color='C'+str(s), s=20, edgecolors='k', linewidth=0.5)


            for stripe, s in zip(stripes, range(1,len(stripes)+1)):

                # mainshock stripe msa fragilities
                current_ax = ax[0][1]
                key = damaged_group + '/' + gm_id + '/' + stripe + '/msa_sa_avg/collapse_fragilities'
                msa_fragilities.loc[stripe,:] = pd.read_hdf(results_filename, key).to_numpy()

                median = msa_fragilities.loc[stripe,'Median']
                beta = msa_fragilities.loc[stripe,'Beta']


                y = np.linspace(0.001, 1, num=100)
                x = stats.lognorm(beta, scale=median).ppf(y)
                label = 'Mainshock: '+ stripe + '\n$IM_{0.5}=$'+'{0:.2f}'.format(median) + ', $\sigma_{ln}=$'+'{0:.2f}'.format(beta)
                _ = current_ax.plot(x, 100*y, label=label, color='C'+str(s), linewidth=fragility_linewidth)

                key = damaged_group + '/' + gm_id + '/' + stripe + '/msa_sa_avg/collapse_matrix'
                msa_stripes = pd.read_hdf(results_filename, key).columns
                collapse_matrix = pd.read_hdf(results_filename, key).to_numpy()
                [n_gms, _] = collapse_matrix.shape
                percent_collapse = 100*np.sum(collapse_matrix, axis=0)/n_gms
                _ = current_ax.scatter(msa_stripes, percent_collapse, color='C'+str(s), edgecolors='k', linewidth=1)

                # mainshock stripe ida fragilities
                current_ax = ax[1][1]
                key = damaged_group + '/' + gm_id + '/' + stripe + '/ida/collapse_fragilities'
                ida_fragilities.loc[stripe,:] = pd.read_hdf(results_filename, key).to_numpy()[1,:] # 1 is for Sa_avg, 0 for Sa(T1)

                median = ida_fragilities.loc[stripe,'Median']
                beta = ida_fragilities.loc[stripe,'Beta']

                y = np.linspace(0.001, 1, num=100)
                x = stats.lognorm(beta, scale=median).ppf(y)
                label = 'Mainshock: '+ stripe + '\n$IM_{0.5}=$'+'{0:.2f}'.format(median) + ', $\sigma_{ln}=$'+'{0:.2f}'.format(beta)
                _ = current_ax.plot(x, 100*y, label=label, color='C'+str(s), linewidth=fragility_linewidth)

                key = damaged_group + '/' + gm_id + '/' + stripe + '/ida/collapse_intensities'
                collapse_intensities = pd.read_hdf(results_filename, key)['Sa_avg'].to_numpy()
                _ = current_ax.scatter(np.sort(collapse_intensities), np.linspace(100/n_gms, 100, num=n_gms, endpoint=True),
                                       color='C'+str(s), s=20, edgecolors='k', linewidth=0.5)


                # post-mainshock idas
                damaged_ida_folder = '/mainshock_damage_results/' + gm_id + '/' + stripe + '/ida/'
                damaged_ida_segments = collect_ida_curves(results_filename, gm_metadata, damaged_ida_folder)
                # adjust to include residual drift
                damaged_edp_folder = '/mainshock_damage_results/' + gm_id + '/' + stripe + '/mainshock_edp/'
                dset = damaged_edp_folder + 'residual_story_drift'
                residual_drift = hf[dset].attrs['max_residual_story_drift']
                damaged_ida_segments[:,0,0] = residual_drift


                if True:
                    current_ax = ax[1][0]
                    ida_plot = LineCollection(damaged_ida_segments, linewidths=1, colors='C'+str(s), linestyle='solid', label=stripe)
                    _ = current_ax.add_collection(ida_plot)
                    _ = current_ax.set_xlim(0,0.1)
                    _ = current_ax.set_ylim(ida_ylim)
                    _ = current_ax.legend(loc='upper left', handlelength=0.5)


            current_ax = ax[0][1]
            legend_title = 'Fragilities via MSA'
            _ = current_ax.legend(title=legend_title, bbox_to_anchor=(1,0.5), loc='center left')
            _ = current_ax.set_ylabel('Probability of Collapse')
            _ = current_ax.set_ylim(0,100)
            _ = current_ax.set_xlabel('Sa$_{avg}$(T$_1$), [g]')
            _ = current_ax.grid('on')

            current_ax = ax[1][1]
            legend_title = 'Fragilities via IDA'
            _ = current_ax.legend(title=legend_title, bbox_to_anchor=(1,0.5), loc='center left')
            _ = current_ax.set_ylabel('Probability of Collapse')
            _ = current_ax.set_ylim(0,100)
            _ = current_ax.set_xlabel('Sa$_{avg}$(T$_1$), [g]')
            _ = current_ax.grid('on')

            current_ax = ax[1][0]
            _ = current_ax.set_xlim(0,0.1)
            _ = current_ax.set_ylim(ida_ylim)
            _ = current_ax.set_ylabel('Sa$_{avg}$(T$_1$), [g]')
            _ = current_ax.set_xlabel('Peak Story Drift Ratio')

            # add reference points to ida
            current_ax = ax[0][0]
            sa_avg_col = gm_metadata.loc[gm_id, 'Intact Collapse Sa_avg']
            stripe_values = sa_avg_col * np.insert(np.array([float(x[:-3]) for x in stripes]), 0, 0)
            x = np.interp(stripe_values, intact_ida_segments[i,:,1], intact_ida_segments[i,:,0])
            color = ['C'+str(k) for k in range(len(stripe_values))]
            _ = current_ax.scatter(x, stripe_values, color=color, zorder=50, s=100)

            fig.tight_layout(rect=[0, 0, 1, 0.95])
            plt.show()


def plot_damaged_ida_per_gm(gm_id, results_filename, gm_metadata, ida_intact_fragility, intact_ida_segments,
                            peak_segments, residual_segments, savefig):
    damaged_group = 'mainshock_damage_results'
    fragility_linewidth = 3

    with h5py.File(results_filename, 'r') as hf:
        i = int(gm_id[2:]) - 1

        ida_fragilities = ida_intact_fragility.copy()
        stripes = list(hf[damaged_group][gm_id].keys())

        if len(stripes) > 1:
            fig, ax = plt.subplots(2, 2, figsize=(15, 10), gridspec_kw={'width_ratios': [0.55, 0.45]})
            suptitle = fig.suptitle('Mainshock: ' + gm_id)
            ax[0][0].get_shared_x_axes().join(ax[0][0], ax[1][0])
            ax[0][0].get_shared_y_axes().join(ax[0][0], ax[1][0])
            ax[0][0].get_shared_y_axes().join(ax[0][0], ax[0][1])

            # intact ida
            current_ax = ax[0][0]
            ida_ylim = [0, 1.1 * np.amax(intact_ida_segments[:, :, 1])]
            full_ida_plot = LineCollection(intact_ida_segments, linewidths=1, colors='lightgray', linestyle='solid')
            _ = current_ax.add_collection(full_ida_plot)
            _ = current_ax.set_xlim(0, 0.1)
            _ = current_ax.set_ylim(ida_ylim)
            _ = current_ax.set_ylabel('Sa$_{avg}$(T$_1$), [g]')
            _ = current_ax.set_xlabel('Peak Story Drift Ratio')

            _ = current_ax.plot(intact_ida_segments[i, :, 0], intact_ida_segments[i, :, 1], color='k', linewidth=2.5,
                                label=gm_id)
            _ = current_ax.legend(loc='upper left')

            if True:
                current_ax = ax[1][0]
                full_ida_plot = LineCollection(intact_ida_segments, linewidths=1, colors='lightgray', linestyle='solid')
                _ = current_ax.add_collection(full_ida_plot)

            # mainshock responses
            current_ax = ax[0][1]
            residual_plot = LineCollection(residual_segments, linewidths=1, colors='dimgray', linestyle='solid')
            _ = current_ax.add_collection(residual_plot)
            peak_plot = LineCollection(peak_segments, linewidths=1, colors='lightgray', linestyle='solid')
            _ = current_ax.add_collection(peak_plot)
            _ = current_ax.plot(residual_segments[i, :, 0], residual_segments[i, :, 1], color='k', linewidth=2.5)
            _ = current_ax.plot(peak_segments[i, :, 0], peak_segments[i, :, 1], color='k', linewidth=2.5)
            legend_elements = [Line2D([0], [0], color='lightgray', label='Peak Drift'),
                               Line2D([0], [0], color='dimgray', label='Residual Drift')]
            _ = current_ax.legend(handles=legend_elements, loc='upper left')
            _ = current_ax.set_xlim(0, 0.06)

            s = 0
            # intact ida fragilities
            current_ax = ax[1][1]
            median = ida_fragilities.loc['Intact', 'Median']
            beta = ida_fragilities.loc['Intact', 'Beta']

            y = np.linspace(0.001, 1, num=100)
            x = stats.lognorm(beta, scale=median).ppf(y)
            label = 'No Mainshock' + '\n$IM_{0.5}=$' + '{0:.2f}'.format(median) + ' $\sigma_{ln}=$' + '{0:.2f}'.format(
                beta)
            _ = current_ax.plot(x, 100 * y, label=label, color='C' + str(s), linewidth=2.5)

            key = '/intact_results/ida/collapse_intensities'
            collapse_intensities = pd.read_hdf(results_filename, key)['Sa_avg'].to_numpy()
            n_gms = len(collapse_intensities)
            _ = current_ax.scatter(np.sort(collapse_intensities),
                                   np.linspace(100 / n_gms, 100, num=n_gms, endpoint=True),
                                   color='C' + str(s), s=20, edgecolors='k', linewidth=0.5)

            for stripe, s in zip(stripes, range(1, len(stripes) + 1)):

                # mainshock stripe ida fragilities
                current_ax = ax[1][1]
                key = damaged_group + '/' + gm_id + '/' + stripe + '/ida/collapse_fragilities'
                ida_fragilities.loc[stripe, :] = pd.read_hdf(results_filename, key).to_numpy()[1,
                                                 :]  # 1 is for Sa_avg, 0 for Sa(T1)

                median = ida_fragilities.loc[stripe, 'Median']
                beta = ida_fragilities.loc[stripe, 'Beta']

                y = np.linspace(0.001, 1, num=100)
                x = stats.lognorm(beta, scale=median).ppf(y)
                label = 'Mainshock: ' + stripe + '\n$IM_{0.5}=$' + '{0:.2f}'.format(
                    median) + ', $\sigma_{ln}=$' + '{0:.2f}'.format(beta)
                _ = current_ax.plot(x, 100 * y, label=label, color='C' + str(s), linewidth=fragility_linewidth)

                key = damaged_group + '/' + gm_id + '/' + stripe + '/ida/collapse_intensities'
                collapse_intensities = pd.read_hdf(results_filename, key)['Sa_avg'].to_numpy()
                n_gms = len(collapse_intensities)
                _ = current_ax.scatter(np.sort(collapse_intensities),
                                       np.linspace(100 / n_gms, 100, num=n_gms, endpoint=True),
                                       color='C' + str(s), s=20, edgecolors='k', linewidth=0.5)

                # post-mainshock idas
                damaged_ida_folder = '/mainshock_damage_results/' + gm_id + '/' + stripe + '/ida/'
                damaged_ida_segments = collect_ida_curves(results_filename, gm_metadata, damaged_ida_folder)
                # adjust to include residual drift
                damaged_edp_folder = '/mainshock_damage_results/' + gm_id + '/' + stripe + '/mainshock_edp/'
                dset = damaged_edp_folder + 'residual_story_drift'
                residual_drift = hf[dset].attrs['max_residual_story_drift']
                damaged_ida_segments[:, 0, 0] = residual_drift

                if True:
                    current_ax = ax[1][0]
                    ida_plot = LineCollection(damaged_ida_segments, linewidths=1, colors='C' + str(s),
                                              linestyle='solid', label=stripe)
                    _ = current_ax.add_collection(ida_plot)
                    _ = current_ax.legend(loc='upper left', handlelength=0.5)

            current_ax = ax[1][1]
            legend_title = 'Fragilities via IDA'
            _ = current_ax.legend(title=legend_title, bbox_to_anchor=(1, 0.5), loc='center left')
            _ = current_ax.set_ylabel('Probability of Collapse')
            _ = current_ax.set_ylim(0, 100)
            _ = current_ax.set_xlabel('Sa$_{avg}$(T$_1$), [g]')
            _ = current_ax.grid('on')

            current_ax = ax[1][0]
            _ = current_ax.set_xlim(0, 0.1)
            _ = current_ax.set_ylim(ida_ylim)
            _ = current_ax.set_ylabel('Sa$_{avg}$(T$_1$), [g]')
            _ = current_ax.set_xlabel('Peak Story Drift Ratio')

            current_ax = ax[0][1]
            _ = current_ax.set_ylabel('Sa$_{avg}$(T$_1$), [g]')
            _ = current_ax.set_xlabel('Story Drift Ratio')

            # add reference points to ida
            current_ax = ax[0][0]
            sa_avg_col = gm_metadata.loc[gm_id, 'Intact Collapse Sa_avg']
            stripe_values = sa_avg_col * np.insert(np.array([float(x[:-3]) for x in stripes]), 0, 0)
            x = np.interp(stripe_values, intact_ida_segments[i, :, 1], intact_ida_segments[i, :, 0])
            color = ['C' + str(k) for k in range(len(stripe_values))]
            _ = current_ax.scatter(x, stripe_values, color=color, zorder=50, s=100)

            # add reference points to mainshocks
            current_ax = ax[0][1]
            y = peak_segments[i, :, 1]
            x = peak_segments[i, :, 0]
            _ = current_ax.scatter(x, y, color=color, zorder=50, s=100)
            x = residual_segments[i, :, 0]
            _ = current_ax.scatter(x, y, color=color, zorder=50, s=100)

            fig.tight_layout(rect=[0, 0, 1, 0.95])
            if savefig:
                figname = gm_id + '_Mainshocks.png'
                plt.savefig(figname, dpi=300)
            plt.show()


def plot_building_at_t(t, edp, columns, beams, plot_scale, ax):
    ax.cla()

    [_, n_pts] = edp.shape
    edp = np.insert(edp, 0, np.zeros((1, n_pts)), axis=0)

    [n_columns, _, _] = columns.shape
    [n_beams, _, _] = beams.shape
    n_stories = int(n_columns - n_beams)
    n_bays = int(n_beams / n_stories)

    columns_t = columns.copy()
    beams_t = beams.copy()

    i_col = 0
    i_beam = 0
    for i_story in range(n_stories):
        for i_end in range(2):
            columns_t[i_col:i_col + n_bays + 2, i_end, 0] = columns[i_col:i_col + n_bays + 2, i_end, 0] + plot_scale * \
                                                                                                          edp[
                                                                                                              i_story + i_end, t]
        i_col = i_col + n_bays + 1

        beams_t[i_beam:i_beam + n_bays + 1, :, 0] = beams[i_beam:i_beam + n_bays + 1, :, 0] + plot_scale * edp[
            i_story + 1, t]
        i_beam = i_beam + n_bays

    column_collection = LineCollection(columns, color='darkgray')
    _ = ax.add_collection(column_collection)

    beam_collection = LineCollection(beams, color='darkgray')
    _ = ax.add_collection(beam_collection)

    column_collection = LineCollection(columns_t, color='k')
    _ = ax.add_collection(column_collection)

    beam_collection = LineCollection(beams_t, color='k')
    _ = ax.add_collection(beam_collection)

    _ = ax.axis('scaled')

    building_height = np.max(columns[:, :, 1])
    building_width = np.max(columns[:, :, 0])
    y_gap = 20
    x_gap = 500
    _ = ax.set_xlim(-x_gap, building_width + x_gap)
    _ = ax.set_ylim(0, building_height + y_gap)
    _ = ax.axis('off')
    # _ = ax.text(building_width / 2, -y_gap, 'Displacement scale: ' + str(plot_scale) + 'x', ha='center', va='top',
    #             fontsize=18)

def plot_hinges(t, edp, joints_x, joints_y, plot_scale, peak_joint_pos, peak_joint_neg, hinge_yield_rotation_positive,
                hinge_cap_rotation_positive, hinge_yield_rotation_negative, hinge_cap_rotation_negative,  ax):
    ## plot hinges takes the rotation of every hinge of a frame around each beam-column joint of coordinates (joint_x by joint_y)
    # in positive and negative direction (peak_joint_pos, peak_joint_neg) and plots the hinges that are in either of
    # the following groups: (yield, cap/3]; (cap/3, cap*2/3], (cap*2/3, cap], or greater than cap
    # INPUTS
    #       t                             = time index to display deformed shape and hinges state
    #       edp                           = 2D displacement array [story, time]
    #       joints_x                      = 2D x-coord of each joint [story, column]
    #       joints_y                      = 2D y-coord of each joint [story, column]
    #       plot_scale                    = scalar to amplify displacement
    #       peak_joint_pos                = 4D array with rotation demand per hinge in each joint [floor, column, hinge_loc, time]
    #                                       hinge_loc is 0: bottom of the joint
    #                                                    1: right
    #                                                    2: top
    #                                                    3: left
    #       peak_joint_neg                = same as before but with maximum negative demand
    #       hinge_yield_rotation_positive = 4D array with yield rotation capacity per hinge in each joint [floor, column, hinge_loc, time]
    #       hinge_cap_rotation_positive
    #       hinge_yield_rotation_negative
    #       hinge_cap_rotation_negative
    #       ax                            = axis to add the hinge plots

    # Retrieve basic info for loops
    n_stories, n_bays = joints_x.shape
    n_stories = n_stories - 1
    n_bays = n_bays - 1

    # Assemble vector with hinges in each state
    disp_t = edp[:, t]  # displacement for deformed shape
    disp_t = np.insert(disp_t, 0, 0, axis=0)  # add the hinge at column base# add zero displacement at base of the column
    dhinge = 4  # plotting delta from joint

    joints_x_yield = np.empty((0, 1))
    joints_y_yield = np.empty((0, 1))

    joints_x_capt = np.empty((0, 1))
    joints_y_capt = np.empty((0, 1))

    joints_x_cap2t = np.empty((0, 1))
    joints_y_cap2t = np.empty((0, 1))

    joints_x_cap = np.empty((0, 1))
    joints_y_cap = np.empty((0, 1))

    for floor_i in range(n_stories + 1):
        disp_curr = disp_t[floor_i] * plot_scale

        for col_i in range(n_bays + 1):

            for hinge_loc in range(4):

                # Read rotation demand of current hinge
                peakPos = peak_joint_pos[floor_i, col_i, hinge_loc, 0]
                peakNeg = -peak_joint_neg[floor_i, col_i, hinge_loc, 0]

                # Read rotation capacity of current hinge
                yieldPos = hinge_yield_rotation_positive[floor_i, col_i, hinge_loc, 0]
                yieldNeg = -hinge_yield_rotation_negative[floor_i, col_i, hinge_loc, 0]
                capPos = hinge_cap_rotation_positive[floor_i, col_i, hinge_loc, 0]
                capNeg = -hinge_cap_rotation_negative[floor_i, col_i, hinge_loc, 0]

                # Plotting position of current hinge
                if hinge_loc == 0:  # Bottom hinge
                    curr_x = joints_x[floor_i, col_i] + disp_curr
                    curr_y = joints_y[floor_i, col_i] - dhinge * plot_scale
                elif hinge_loc == 1:  # Right hinge
                    curr_x = joints_x[floor_i, col_i] + disp_curr + dhinge * plot_scale
                    curr_y = joints_y[floor_i, col_i]
                elif hinge_loc == 2:  # Top hinge
                    curr_x = joints_x[floor_i, col_i] + disp_curr
                    curr_y = joints_y[floor_i, col_i] + dhinge * plot_scale
                else:  # Left hinge
                    curr_x = joints_x[floor_i, col_i] + disp_curr - dhinge * plot_scale
                    curr_y = joints_y[floor_i, col_i]

                if yieldPos != 0:
                    if (peakPos > yieldPos and peakPos <= capPos / 3) or (peakNeg > yieldNeg and peakNeg <= capNeg / 3):
                        # Between yield and capping/3
                        joints_x_yield = np.append(joints_x_yield, curr_x)
                        joints_y_yield = np.append(joints_y_yield, curr_y)

                    elif (peakPos > capPos / 3 and peakPos <= capPos * 2 / 3) or (
                            peakNeg > capNeg / 3 and peakNeg <= capNeg * 2 / 3):
                        # Between capping/3 and capping*2/3
                        joints_x_capt = np.append(joints_x_capt, curr_x)
                        joints_y_capt = np.append(joints_y_capt, curr_y)

                    elif (peakPos > capPos * 2 / 3 and peakPos <= capPos) or (
                            peakNeg > capNeg * 2 / 3 and peakNeg <= capNeg):
                        # Between capping*2/3 and capping
                        #                         print('yieldPos='+str(yieldPos)+';capPos='+str(capPos))
                        #                         print('floor='+str(floor_i)+';column='+str(col_i)+';hinge='+str(hinge_loc))
                        joints_x_cap2t = np.append(joints_x_cap2t, curr_x)
                        joints_y_cap2t = np.append(joints_y_cap2t, curr_y)
                    elif (peakPos > capPos) or (peakNeg > capNeg):
                        # beyond capping point
                        #                         print('yieldPos='+str(yieldPos)+';capPos='+str(capPos))
                        #                         print('floor='+str(floor_i)+';column='+str(col_i)+';hinge='+str(hinge_loc))
                        joints_x_cap = np.append(joints_x_cap, curr_x)
                        joints_y_cap = np.append(joints_y_cap, curr_y)

    # Plot hinges with appropiate colors
    _ = ax.plot(joints_x_yield, joints_y_yield, 'o', color='m')
    _ = ax.plot(joints_x_capt, joints_y_capt, 'o', color='c')
    _ = ax.plot(joints_x_cap2t, joints_y_cap2t, 'o', color='r')
    _ = ax.plot(joints_x_cap, joints_y_cap, 'o', color='b')

def plot_mainshock_damage_visual(displacement, periods, spectrum, acc, dt, n_pts, time_series, story_idx, story_displacement,
                                 story_drift, column_geometry, beam_geometry, joints_x, joints_y,peak_joint_pos,
                                 peak_joint_neg, hinge_yield_rotation_positive, hinge_cap_rotation_positive,
                                 hinge_yield_rotation_negative, hinge_cap_rotation_negative):

    ax = list()

    fig = plt.figure(figsize=(15, 10))
    ax.append(plt.subplot2grid((3, 3), (0, 0), rowspan=1, colspan=2))
    ax.append(plt.subplot2grid((3, 3), (1, 0), rowspan=1, colspan=2))
    ax.append(plt.subplot2grid((3, 3), (2, 0), rowspan=1, colspan=2))
    ax.append(plt.subplot2grid((3, 3), (0, 2), rowspan=2, colspan=1))
    ax.append(plt.subplot2grid((3, 3), (2, 2), rowspan=1, colspan=1))

    for i in [1, 2]:
        ax[0].get_shared_x_axes().join(ax[0], ax[i])

    color = 'tab:blue'

    current_ax = ax[4]
    _ = current_ax.plot(periods, spectrum, color=color)
    _ = current_ax.set_xlim(0, 5)
    _ = current_ax.set_ylabel('Spectral Acceleration [g]')
    _ = current_ax.set_xlabel('Period')

    current_ax = ax[2]
    _ = current_ax.plot(np.arange(n_pts) * dt, acc, color=color)
    _ = current_ax.set_xlim(left=0)
    ylim = 1.1 * np.max(np.abs(acc))
    _ = current_ax.set_ylim((-ylim, ylim))
    _ = current_ax.set_ylabel('Ground Acceleration [g]')
    _ = current_ax.set_xlabel('Seconds')

    edp_list = [story_displacement, story_drift]
    edp_name = ['Displacement [in]', 'Story Drift Ratio [%]']
    level_name = [str(story_idx + 1) + 'th Floor', str(story_idx) + 'th Story']

    for edp, i in zip(edp_list, range(len(edp_list))):
        current_ax = ax[i]
        _ = current_ax.plot(time_series, edp, color=color)
        ylim = 1.1 * np.max(np.abs(edp))
        _ = current_ax.set_ylim((-ylim, ylim))
        ylabel = level_name[i] + '\n' + edp_name[i]
        _ = current_ax.set_ylabel(ylabel)

    for i in range(5):
        if i != 3:
            _ = ax[i].grid('on')

    current_ax = ax[0]
    _ = current_ax.set_xlim((0, time_series[-1]))

    current_ax = ax[3]
    t = len(time_series) - 1
    plot_scale = 10
    plot_building_at_t(t, displacement, column_geometry, beam_geometry, plot_scale, current_ax)
    plot_hinges(t, displacement, joints_x, joints_y, plot_scale, peak_joint_pos, peak_joint_neg,
                hinge_yield_rotation_positive, hinge_cap_rotation_positive, hinge_yield_rotation_negative,
                hinge_cap_rotation_negative, current_ax)
    plt.tight_layout()
    plt.show()


def plot_increasing_two_bin_threshold(edp_type, cutoff_list, medians, betas, n_instances, intact_median, intact_beta):
    fig = plt.figure(figsize=(12, 6))

    intact_color = 'k'
    above_color = 'tab:red'
    below_color = 'tab:blue'
    s = 100

    ax = list()
    n_col = 3
    n_row = 100
    col_scale = 0.4
    col_width = int(n_row * col_scale)
    col_start = int((0.5 - col_scale) / 2)
    col1 = col_start
    col2 = int(n_row * 0.5 + col_start)
    ax.append(plt.subplot2grid((n_col, n_row), (0, col2), rowspan=1, colspan=col_width))
    ax.append(plt.subplot2grid((n_col, n_row), (1, col1), rowspan=n_col - 1, colspan=col_width))
    ax.append(plt.subplot2grid((n_col, n_row), (1, col2), rowspan=n_col - 1, colspan=col_width))

    for i in [1, 2]:
        ax[0].get_shared_x_axes().join(ax[0], ax[i])

    cutoff_list = cutoff_list * 100

    xlabel_prefix = 'Damage Indicator Threshold:\n'
    if edp_type == 'peak':
        xlabel = xlabel_prefix + 'Peak Story Drift Ratio [%]'
    elif edp_type == 'residual':
        xlabel = xlabel_prefix + 'Residual Story Drift Ratio [%]'

    current_ax = ax[0]
    total_instances = np.sum(n_instances[0, :])
    bin_percents = 100 * n_instances[:, 0] / total_instances

    _ = current_ax.fill_between(cutoff_list, 100, bin_percents, color=above_color, alpha=0.5)
    _ = current_ax.fill_between(cutoff_list, bin_percents, 0, color=below_color, alpha=0.5)
    _ = current_ax.set_xlim(0, max(cutoff_list))
    _ = current_ax.set_ylim(0, 100)
    _ = current_ax.set_ylabel('Above and\nBelow [%]')
    _ = current_ax.set_xticks([])
    _ = current_ax.set_yticks([0, 100])

    current_ax = ax[1]
    _ = current_ax.scatter(0, intact_median, s=s, color=intact_color, label='Intact Fragility', zorder=10)
    _ = current_ax.plot(cutoff_list, medians[:, 1], color=above_color, label='Above DI Threshold')
    _ = current_ax.plot(cutoff_list, medians[:, 0], color=below_color, label='Below DI Threshold')
    _ = current_ax.set_ylabel('Median Collapse Capacity')
    _ = current_ax.set_ylim(bottom=0)
    _ = current_ax.set_xlabel(xlabel)
    _ = current_ax.legend(bbox_to_anchor=(0.5, 1.1), loc='lower center')

    current_ax = ax[2]
    _ = current_ax.scatter(0, intact_beta, s=s, color=intact_color, zorder=10)
    _ = current_ax.plot(cutoff_list, betas[:, 1], color=above_color)
    _ = current_ax.plot(cutoff_list, betas[:, 0], color=below_color)
    _ = current_ax.set_ylabel('Fragility Dispersion')
    _ = current_ax.set_xlabel(xlabel)


def plot_two_bin_fragilities(intact_median, medians, intact_beta, betas, bin_max, threshold):
    i = np.where(bin_max[:, 0] == threshold)

    intact_color = 'k'
    above_color = 'tab:red'
    below_color = 'tab:blue'

    colors = [intact_color, below_color, above_color]
    label_prefixes = ['Intact', 'DI <= ' + '{0:.2f}'.format(100 * threshold) + '%',
                      'DI > ' + '{0:.2f}'.format(100 * threshold) + '%']

    median = np.append(intact_median, medians[i, :])
    beta = np.append(intact_beta, betas[i, :])
    threshold = bin_max[i, 0]

    fig, current_ax = plt.subplots(1, 1)
    _ = current_ax.plot()

    for j in range(3):
        y = np.linspace(0.001, 1, num=100)
        x = stats.lognorm(beta[j], scale=median[j]).ppf(y)
        label = label_prefixes[j] + '\n$IM_{0.5}=$' + '{0:.2f}'.format(
            median[j]) + ', $\sigma_{ln}=$' + '{0:.2f}'.format(beta[j])
        _ = current_ax.plot(x, 100 * y, label=label, color=colors[j])

    _ = plt.legend(bbox_to_anchor=(1, 0.5), loc='center left')
    _ = current_ax.set_xlabel('Sa$_{avg}$(T$_1$), [g]')
    _ = current_ax.set_xlim(left=0)
    _ = current_ax.set_ylabel('Probability of Collapse')
    _ = current_ax.set_ylim(0, 100)

    _ = plt.show()


def plot_multi_bin_fragilities(intact_fragility, fragilities, edp_k_raghunandan, edp_k_burton, p_values, edp_type):

    intact_median = intact_fragility['Median'].values

    max_edps = 100 * fragilities['Max EDP'].values
    fragilities = fragilities.loc[:, ['Median', 'Beta']]
    fragilities = pd.concat([intact_fragility, fragilities])

    target_edp = 100 * fragilities.index.values
    target_edp[-1] = max_edps[-1]

    medians = fragilities['Median'].values
    betas = fragilities['Beta'].values
    ks = medians / intact_median

    if edp_type == 'peak':
        edp_label = 'Peak Drift '
    elif edp_type == 'residual':
        edp_label = 'Res. Drift '

    fig, ax = plt.subplots(1, 2, figsize=(15, 5))

    i = 0
    for edp, s in zip(fragilities.index, range(len(fragilities))):

        median = fragilities.loc[edp, 'Median']
        beta = fragilities.loc[edp, 'Beta']

        y = np.linspace(0.001, 1, num=100)
        x = stats.lognorm(beta, scale=median).ppf(y)

        if s == 0:
            label = 'Intact' + '\n$IM_{0.5}=$' + '{0:.2f}'.format(median) + ' $\sigma_{ln}=$' + '{0:.2f}'.format(beta)
        elif s != len(fragilities) - 1:
            label = edp_label + '<= ' + str(100 * edp) + '%' + '\n$IM_{0.5}=$' + '{0:.2f}'.format(
                median) + ', $\sigma_{ln}=$' + '{0:.2f}'.format(beta)
        else:
            label = edp_label + '> ' + str(100 * fragilities.index[s - 1]) + '%' + '\n$IM_{0.5}=$' + '{0:.2f}'.format(
                median) + ', $\sigma_{ln}=$' + '{0:.2f}'.format(beta)
        _ = ax[i].plot(x, 100 * y, label=label, color='C' + str(s), linewidth=2.5)

    _ = ax[i].legend(bbox_to_anchor=(1, 0.5), loc='center left')
    _ = ax[i].set_ylabel('Probability of Collapse')
    _ = ax[i].set_ylim(0, 100)
    _ = ax[i].set_xlabel('Sa$_{avg}$(T$_1$), [g]')
    _ = ax[i].grid('on')

    i = 1
    for edp, s in zip(fragilities.index, range(len(fragilities))):

        median = fragilities.loc[edp, 'Median']
        beta = fragilities.loc[edp, 'Beta']

        y = np.linspace(0.001, 1, num=100)
        x = stats.lognorm(beta, scale=median).ppf(y) / intact_median
        if s == 0:
            label = 'Intact' + '\n$IM_{0.5}=$' + '{0:.2f}'.format(median) + ' $\sigma_{ln}=$' + '{0:.2f}'.format(beta)
        elif s != len(fragilities) - 1:
            label = edp_label + '<= ' + str(100 * edp) + '%' + '\n$IM_{0.5}=$' + '{0:.2f}'.format(
                median) + ', $\sigma_{ln}=$' + '{0:.2f}'.format(beta)
        else:
            label = edp_label + '> ' + str(100 * fragilities.index[s - 1]) + '%' + '\n$IM_{0.5}=$' + '{0:.2f}'.format(
                median) + ', $\sigma_{ln}=$' + '{0:.2f}'.format(beta)
        _ = ax[i].plot(x, 100 * y, label=label, color='C' + str(s), linewidth=2.5)

        _ = ax[i].scatter(median / intact_median, 50, color='C' + str(s), s=50, edgecolor='k', zorder=10)
        for p in p_values:
            p_x = (stats.lognorm(beta, scale=median).ppf(p)) / intact_median
            _ = ax[i].scatter(p_x, 100 * p, color='k', marker='|', s=100, zorder=10)

    _ = ax[i].set_ylabel('Probability of Collapse')
    _ = ax[i].set_ylim(0, 100)
    _ = ax[i].set_xlabel('Sa$_{avg}$(T$_1$) / $\widehat{Sa}_{avg}$(T$_1)_{intact}$')
    _ = ax[i].grid('on')

    plt.tight_layout()

    fig, ax = plt.subplots(1, 2, figsize=(15, 5))

    color = 'lightgray'
    alpha = 1
    _ = ax[0].scatter(edp_k_raghunandan['EDP'], edp_k_raghunandan['kappa'], facecolor='none', edgecolor=color,
                      alpha=alpha, zorder=-10)
    _ = ax[1].scatter(edp_k_burton['EDP'], edp_k_burton['kappa'], facecolor='none', edgecolor=color, alpha=alpha,
                      zorder=-10)

    for i in range(2):
        colors = ['C' + str(j) for j in range(len(target_edp))]
        _ = ax[i].scatter(target_edp, ks, color=colors, edgecolor='k', zorder=50, s=150)

        upper_bound = [stats.lognorm(betas[i], scale=medians[i]).ppf(p_values[0]) for i in
                       range(len(medians))] / intact_median - ks
        lower_bound = ks - [stats.lognorm(betas[i], scale=medians[i]).ppf(p_values[1]) for i in
                            range(len(medians))] / intact_median

        _ = ax[i].errorbar(target_edp, ks, yerr=[lower_bound, upper_bound], fmt='o', color='k', capsize=4, capthick=2,
                           zorder=40)

        _ = ax[i].plot([max_edps] * 2, [0, 1.2], color='dimgray', linestyle='--', zorder=0)

        if edp_type == 'peak':
            _ = ax[i].set_xlabel('Peak Story Drift Ratio [%]')
        elif edp_type == 'residual':
            _ = ax[i].set_xlabel('Residual Story Drift Ratio [%]')

        if i == 0:
            ylabel = 'Reduction in Collapse Capacity,\n' + r'$\kappa=\frac{Sa_{col, damaged}}{Sa_{col, intact}}$'
            _ = ax[i].set_ylabel(ylabel)
            _ = ax[i].set_title('Comparison with Raghunandan et al. 2015')
        if i == 1:
            ylabel = 'Reduction in Collapse Capacity,\n' + r'$\kappa=\frac{\widehat{Sa}_{col, damaged}}{\widehat{Sa}_{col, intact}}$'
            _ = ax[i].set_ylabel(ylabel)
            _ = ax[i].set_title('Comparison with Burton et al. 2018')

    plt.tight_layout()
    _ = plt.show()
