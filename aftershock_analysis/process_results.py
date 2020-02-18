from .base import *
from .collect_nrha_results import *


def k_by_damage_instance(results_filename, edp_results, fragility_type, ax):  # add fragility type (msa or ida) as an input

    gm_ids = edp_results.index
    scales = edp_results.columns

    n_gms = len(gm_ids)
    n_scales = len(scales)

    edp = np.zeros((n_gms * n_scales))
    k = np.zeros_like(edp)

    if fragility_type == 'msa_sa_avg':
        fragility_key = '/msa_sa_avg/collapse_fragilities'
        fragility_idx = 0
    elif fragility_type == 'ida_sa_avg':
        fragility_key = '/ida/collapse_fragilities'
        fragility_idx = 1
    key = 'intact_results' + fragility_key
    intact_median = pd.read_hdf(results_filename, key).to_numpy()[fragility_idx][0]

    i = 0
    for gm_id in gm_ids:
        for scale in scales:
            edp[i] = 100 * edp_results.loc[gm_id, scale]

            key = 'mainshock_damage_results/' + gm_id + '/' + str(scale) + 'Col' + fragility_key
            median = pd.read_hdf(results_filename, key).to_numpy()[fragility_idx][0]
            k[i] = median / intact_median

            i = i + 1

    edp_k = pd.DataFrame({'EDP': edp, 'kappa': k})

    if ax is not None:
        _ = ax.scatter(edp, k, facecolor='none', edgecolor='tab:blue')
        _ = ax.set_xlabel('Max Interstory Drift Ratio [%]')
        ylabel = 'Reduction in Collapse Capacity,\n' + r'$\kappa=\frac{\widehat{Sa}_{col, damaged}}{\widehat{Sa}_{col, intact}}$'
        _ = ax.set_ylabel(ylabel)
        _ = ax.set_title('Comparison with Burton et al. 2018')

    return edp_k


def k_by_damage_instance_and_gm(results_filename, edp_results, ax):
    gm_ids = edp_results.index
    scales = edp_results.columns

    n_gms = len(gm_ids)
    n_scales = len(scales)

    edp = np.zeros((n_gms * n_scales * n_gms))
    k = np.zeros_like(edp)

    key = 'intact_results/ida/collapse_intensities'
    intact_collapse_intensities = pd.read_hdf(results_filename, key)['Scale Factor']

    i = 0
    for gm_id in gm_ids:
        for scale in scales:
            start_i = i * n_gms
            end_i = start_i + n_gms
            edp[start_i:end_i] = 100 * edp_results.loc[gm_id, scale]

            if end_i == 7784:
                print(edp.shape)

            key = 'mainshock_damage_results/' + gm_id + '/' + str(scale) + 'Col' + '/ida/collapse_intensities'
            collapse_intensities = pd.read_hdf(results_filename, key)['Scale Factor']
            k[start_i:end_i] = collapse_intensities / intact_collapse_intensities

            i = i + 1

    edp_k = pd.DataFrame({'EDP': edp, 'kappa': k})

    if ax is not None:
        _ = ax.scatter(edp, k, facecolor='none', edgecolor='tab:blue')
        _ = ax.set_xlabel('Max Interstory Drift Ratio [%]')
        ylabel = 'Reduction in Collapse Capacity,\n' + r'$\kappa=\frac{Sa_{col, damaged}}{Sa_{col, intact}}$'
        _ = ax.set_ylabel(ylabel)
        _ = ax.set_title('Comparison with Raghunandan et al. 2015')

    return edp_k


def fragility_by_edp(results_filename, edp_results, edp_categories, edp_cutoffs, fragility_type):
    gm_ids = edp_results.index
    scales = edp_results.columns
    df = pd.DataFrame(index=gm_ids, columns=scales)

    edp_results = edp_results.to_numpy()

    n_categories = len(edp_categories)

    if fragility_type == 'msa_sa_avg':
        gm_scale_group = 'mainshock_damage_results/' + gm_ids[0] + '/' + str(scales[0]) + 'Col'
        key = gm_scale_group + '/msa_sa_avg/collapse_matrix'
        stripes = pd.read_hdf(results_filename, key=key).columns

        fragilities = pd.DataFrame(columns=['Median', 'Beta', 'Min EDP', 'Max EDP', 'N Damaged Instances'],
                                   dtype='float64')
        for i in range(n_categories):
            [gm_idx, scale_idx] = np.where((edp_results > edp_cutoffs[i]) & (edp_results <= edp_cutoffs[i + 1]))
            df.iloc[gm_idx, scale_idx] = edp_categories[i]
            n_instances = len(gm_idx)

            total_collapse_matrix = pd.DataFrame(columns=stripes, dtype='float64')
            for i_gm, i_scale in zip(gm_idx, scale_idx):
                gm_id = gm_ids[i_gm]
                scale_name = str(scales[i_scale]) + 'Col'

                gm_scale_group = 'mainshock_damage_results/' + gm_id + '/' + scale_name
                key = gm_scale_group + '/msa_sa_avg/collapse_matrix'
                collapse_matrix = pd.read_hdf(results_filename, key=key)
                total_collapse_matrix = pd.concat([total_collapse_matrix, collapse_matrix], axis=0)

            if n_instances > 0:
                fragilities.loc[edp_categories[i], ['Median', 'Beta']] = compute_msa_fragility(total_collapse_matrix,
                                                                                           plot=False)
            fragilities.loc[edp_categories[i], 'Min EDP'] = edp_cutoffs[i]
            fragilities.loc[edp_categories[i], 'Max EDP'] = edp_cutoffs[i + 1]
            fragilities.loc[edp_categories[i], 'N Damaged Instances'] = n_instances

    elif fragility_type == 'ida_sa_avg':
        gm_scale_group = 'mainshock_damage_results/' + gm_ids[0] + '/' + str(scales[0]) + 'Col'
        key = gm_scale_group + '/ida/collapse_intensities'
        columns = pd.read_hdf(results_filename, key=key).columns

        fragilities = pd.DataFrame(columns=['Median', 'Beta', 'Max EDP', 'Min EDP', 'N Damaged Instances'],
                                   dtype='float64')
        for i in range(n_categories):
            [gm_idx, scale_idx] = np.where((edp_results > edp_cutoffs[i]) & (edp_results <= edp_cutoffs[i + 1]))
            df.iloc[gm_idx, scale_idx] = edp_categories[i]
            n_instances = len(gm_idx)

            total_collapse_intensities = pd.DataFrame(columns=columns, dtype='float64')
            for i_gm, i_scale in zip(gm_idx, scale_idx):
                gm_id = gm_ids[i_gm]
                scale_name = str(scales[i_scale]) + 'Col'

                gm_scale_group = 'mainshock_damage_results/' + gm_id + '/' + scale_name
                key = gm_scale_group + '/ida/collapse_intensities'
                collapse_intensities = pd.read_hdf(results_filename, key=key)
                total_collapse_intensities = pd.concat([total_collapse_intensities, collapse_intensities], axis=0)

            im = 'Sa_avg'
            if n_instances > 0:
                fragilities.loc[edp_categories[i], ['Median', 'Beta']] = compute_ida_fragility(
                    total_collapse_intensities[im], plot=False)
            fragilities.loc[edp_categories[i], 'Min EDP'] = edp_cutoffs[i]
            fragilities.loc[edp_categories[i], 'Max EDP'] = edp_cutoffs[i + 1]
            fragilities.loc[edp_categories[i], 'N Damaged Instances'] = n_instances

    return fragilities
