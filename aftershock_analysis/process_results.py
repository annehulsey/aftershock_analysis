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

        ida_curve = ida_curve.loc[:, ['Sa_avg', 'Interstory Drift Ratio (max)']].to_numpy()
        collapse = ida_curve[-1, 0]
        ida_segments[i, :, 1] = collapse

        ida_segments[i, 0:len(ida_curve), 1] = ida_curve[:, 0]
        ida_segments[i, 0:len(ida_curve), 0] = ida_curve[:, 1]

    return ida_segments


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
            _ = current_ax.set_xlabel('Maximum Interstory Drift Ratio')

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
            _ = current_ax.set_xlabel('Maximum Interstory Drift Ratio')

            # add reference points to ida
            current_ax = ax[0][0]
            sa_avg_col = gm_metadata.loc[gm_id, 'Intact Collapse Sa_avg']
            stripe_values = sa_avg_col * np.insert(np.array([float(x[:-3]) for x in stripes]), 0, 0)
            x = np.interp(stripe_values, intact_ida_segments[i,:,1], intact_ida_segments[i,:,0])
            color = ['C'+str(k) for k in range(len(stripe_values))]
            _ = current_ax.scatter(x, stripe_values, color=color, zorder=50, s=100)

            fig.tight_layout(rect=[0, 0, 1, 0.95])
            plt.show()