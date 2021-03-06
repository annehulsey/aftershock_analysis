from .base import *


def compute_di_deterministic(isBeam, peakPos, peakNeg, yieldPos, capPos, yieldNeg, capNeg):
    # Function to compute hinge damage index (di)
    if isBeam:
        # is a beam
        if (peakPos > yieldPos and peakPos <= capPos / 3) or (peakNeg > yieldNeg and peakNeg <= capNeg / 3):
            # Between yield and capping/3
            di = 1  # 0.5

        elif (peakPos > capPos / 3 and peakPos <= capPos * 2 / 3) or (
                peakNeg > capNeg / 3 and peakNeg <= capNeg * 2 / 3):
            # Between capping/3 and capping*2/3
            di = 2  # 1.0

        elif (peakPos > capPos * 2 / 3 and peakPos <= capPos) or (
                peakNeg > capNeg * 2 / 3 and peakNeg <= capNeg):
            # Between capping*2/3 and capping
            di = 3.0  # 2.0

        elif (peakPos > capPos) or (peakNeg > capNeg):
            # beyond capping point
            di = 4.0
        else:
            # no damage
            di = 0

    else:
        # is a column
        if (peakPos > yieldPos and peakPos <= capPos / 3) or (peakNeg > yieldNeg and peakNeg <= capNeg / 3):
            # Between yield and capping/3
            di = 1.0  # 2.0

        elif (peakPos > capPos / 3 and peakPos <= capPos * 2 / 3) or (
                peakNeg > capNeg / 3 and peakNeg <= capNeg * 2 / 3):
            # Between capping/3 and capping*2/3
            di = 2  # 2.5

        elif (peakPos > capPos * 2 / 3 and peakPos <= capPos) or (
                peakNeg > capNeg * 2 / 3 and peakNeg <= capNeg):
            # Between capping*2/3 and capping
            di = 3.0  # 3.0

        elif (peakPos > capPos) or (peakNeg > capNeg):
            # beyond capping point
            di = 4.0  # 4.0
        else:
            # no damage
            di = 0

    return di


def get_FEMAP58_fragility(isBeam, rot_cap):
    # Build fragilities based on FEMA P58 method using EDPs
    #
    # INPUTS
    #   isBeam  = bool to use different definition of fragility medians based on capping rotation
    #             for beams and columns
    #   rot_cap = capping rotation of the element
    #
    # OUTPUTS
    #   fragilities = panda DataFrame with medians and betas of the damage fragility for the element

    if isBeam:
        fragilities = pd.DataFrame({'median': np.array([0.3, 0.7, 1.0]) * rot_cap,
                                    'beta': np.array([1, 1, 1]) * 0.4}, index=[1, 2, 3])
    else:
        fragilities = pd.DataFrame({'median': np.array([0.25, 0.55, 0.8]) * rot_cap,
                                    'beta': np.array([1, 1, 1]) * 0.4}, index=[1, 2, 3])

    return fragilities


def compute_di_prob(peakPos, peakNeg, pos_fragilities, neg_fragilities):
    # Function to compute hinge damage index (di) accounting for fragilities but assuming independende of damage
    #
    # INPUTS
    #    peakPos          = float edp demand in positive direction
    #    peakNeg          = float edp demand in negative direction
    #    pos_fragilitieS  = DataFrame with damage fragilities for positive demand
    #    neg_fragilities  = DataFrame with damage fragilities for negative demand
    #
    # OUTPUTS
    #    di               = Expected value of the di
    #

    n_damage_states = len(pos_fragilities['median'])

    cumm_prob = np.zeros(n_damage_states)
    weights = np.zeros(n_damage_states + 1)
    weights_neg = np.zeros(n_damage_states + 1)

    # Compute di for positive peak demand
    for ds_i in range(n_damage_states):
        median = pos_fragilities['median'][ds_i + 1]
        sigma_ln = pos_fragilities['beta'][ds_i + 1]
        cumm_prob[ds_i] = stats.lognorm.cdf(peakPos, sigma_ln, 0, median)

    for ds_i in range(n_damage_states + 1):
        if ds_i == 0:
            weights[ds_i] = 1 - cumm_prob[ds_i]
        elif ds_i <= n_damage_states - 1:
            weights[ds_i] = cumm_prob[ds_i - 1] - cumm_prob[ds_i]
        else:
            weights[ds_i] = cumm_prob[ds_i - 1]

    di = np.dot([0, 1, 2, 3], weights)

    # Compute di for negative peak demand
    for ds_i in range(n_damage_states):
        median = neg_fragilities['median'][ds_i + 1]
        sigma_ln = neg_fragilities['beta'][ds_i + 1]
        cumm_prob[ds_i] = stats.lognorm.cdf(peakNeg, sigma_ln, 0, median)

    for ds_i in range(n_damage_states + 1):
        if ds_i == 0:
            weights_neg[ds_i] = 1 - cumm_prob[ds_i]
        elif ds_i <= 2:
            weights_neg[ds_i] = cumm_prob[ds_i - 1] - cumm_prob[ds_i]
        else:
            weights_neg[ds_i] = cumm_prob[ds_i - 1]

    di_neg = np.dot([0, 1, 2, 3], weights_neg)

    if di_neg > di:
        di = di_neg
        weights = weights_neg

    return di, weights


def compute_FDI(floor_i, peak_joint_pos, peak_joint_neg, hinge_yield_rotation_positive, hinge_cap_rotation_positive,
                hinge_yield_rotation_negative, hinge_cap_rotation_negative, inspector_type):
    # Function to compute floor damage index (FDI)
    # INPUTS
    #    floor_i                            = Index (between 0 and n_stories+1) of the current floor
    #    peak_joint_pos                     = 4D array with the positive peak rotation in each hinge [floor_i,col_i,hinge_loc, 0]
    #                                         hinge_loc = 0 bottom; 1: right; 2: top; 3: left of joint at floor_i, col_i
    #    peak_joint_neg                     = 4D array with the negative peak rotation in each hinge [floor_i,col_i,hinge_loc, 0]
    #                                         hinge_loc = 0 bottom; 1: right; 2: top; 3: left of joint at floor_i, col_i
    #    hinge_yield_rotation_positive      = 4D array with the pos yield rotation capacity in each hinge [floor_i,col_i,hinge_loc, 0]
    #    hinge_cap_rotation_positive        = 4D array with the pos rotation at capping in each hinge [floor_i,col_i,hinge_loc, 0]
    #    hinge_yield_rotation_negative      = 4D array with the neg yield rotation capacity in each hinge [floor_i,col_i,hinge_loc, 0]
    #    hinge_cap_rotation_negative        = 4D array with the neg rotation at capping in each hinge [floor_i,col_i,hinge_loc, 0]
    #    inspector_type                     = 'Deterministic' to use deterministic fragilities to assign damage state
    #                                       = 'Probablistic' to use lognormal fragilities to assign damage state
    # OUTPUTS
    #    FDI       = Floor Damage Index for floor_i [0.0 <= FDI <= 1.0 all beams and columns in last damage state]
    #    FDIbeams  = Floor Damage Index accounting for beams only [0.0 <= FDIbeams <= 1.0 all beams in last damage state]
    #    FDIcol    = Floor Damage Index accounting for columns only [0.0 <= FDIcol <= 1.0 all columns in last damage state]
    #

    # Basic geometrical data
    n_stories, n_bays, _, _ = hinge_yield_rotation_positive.shape
    n_stories = n_stories - 1
    n_bays = n_bays - 1

    # Number of potential hinges at given story
    beam_hinges_per_story = n_bays * 2
    if floor_i == 0 or floor_i == n_stories + 1:
        # at first and last floors only columns on one side
        column_hinges_per_story = (n_bays + 1)
    else:
        # al all other floors columns on top and bottom
        column_hinges_per_story = (n_bays + 1) * 2

    # Initialize array to store di of hinges in the floor
    floor_beam_di = np.zeros([beam_hinges_per_story])
    floor_column_di = np.zeros([column_hinges_per_story])

    # index for beam and columns arrays
    beam_hinge_i = 0
    column_hinge_i = 0

    # Loop on all hinges
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

            # hinge damage index (di)
            if yieldPos != 0:  # and peakPos != 0 and peakNeg != 0:
                if hinge_loc == 1 or hinge_loc == 3:
                    isBeam = True
                    if inspector_type == 'Deterministic':
                        floor_beam_di[beam_hinge_i] = compute_di_deterministic(isBeam, peakPos, peakNeg, yieldPos,
                                                                               capPos,
                                                                               yieldNeg, capNeg)
                    else:
                        pos_fragilities = get_FEMAP58_fragility(isBeam, capPos)
                        neg_fragilities = get_FEMAP58_fragility(isBeam, capNeg)
                        floor_beam_di[beam_hinge_i], _ = compute_di_prob(peakPos, peakNeg, pos_fragilities,
                                                                         neg_fragilities)

                    beam_hinge_i = beam_hinge_i + 1
                else:
                    isBeam = False
                    if inspector_type == 'Deterministic':
                        floor_column_di[column_hinge_i] = compute_di_deterministic(isBeam, peakPos, peakNeg, yieldPos,
                                                                                   capPos, yieldNeg, capNeg)
                    else:
                        pos_fragilities = get_FEMAP58_fragility(isBeam, capPos)
                        neg_fragilities = get_FEMAP58_fragility(isBeam, capNeg)
                        floor_column_di[column_hinge_i], _ = compute_di_prob(peakPos, peakNeg, pos_fragilities,
                                                                             neg_fragilities)

                    column_hinge_i = column_hinge_i + 1

    if inspector_type == 'Deterministic':
        n_ds = 4
    else:
        n_ds = 3

    FDI = (0.5 * sum(floor_beam_di) / n_ds / beam_hinges_per_story + 0.5 * sum(
        floor_column_di) / n_ds / column_hinges_per_story)
    FDIbeams = sum(floor_beam_di) / n_ds / beam_hinges_per_story
    FDIcol = sum(floor_column_di) / n_ds / column_hinges_per_story

    return FDI, FDIbeams, FDIcol

def compute_building_DI(peak_joint_pos, peak_joint_neg, hinge_yield_rotation_positive, hinge_cap_rotation_positive,
                        hinge_yield_rotation_negative, hinge_cap_rotation_negative, inspector_type):
    # compute_building_DI computes a Damage Indicator that captures the level of damage of all the components (beams
    # and columns of a frame). The damage state of each component is inferred using deterministic or probabilistic
    # component fragilities in all cases assuming that damage states for each component are mutually independent
    #
    # INPUTS
    #    peak_joint_pos                     = 4D array with the positive peak rotation in each hinge [floor_i,col_i,hinge_loc, 0]
    #                                         hinge_loc = 0 bottom; 1: right; 2: top; 3: left of joint at floor_i, col_i
    #    peak_joint_neg                     = 4D array with the negative peak rotation in each hinge [floor_i,col_i,hinge_loc, 0]
    #                                         hinge_loc = 0 bottom; 1: right; 2: top; 3: left of joint at floor_i, col_i
    #    hinge_yield_rotation_positive      = 4D array with the pos yield rotation capacity in each hinge [floor_i,col_i,hinge_loc, 0]
    #    hinge_cap_rotation_positive        = 4D array with the pos rotation at capping in each hinge [floor_i,col_i,hinge_loc, 0]
    #    hinge_yield_rotation_negative      = 4D array with the neg yield rotation capacity in each hinge [floor_i,col_i,hinge_loc, 0]
    #    hinge_cap_rotation_negative        = 4D array with the neg rotation at capping in each hinge [floor_i,col_i,hinge_loc, 0]
    #    inspector_type                     = 'Deterministic' to use deterministic fragilities to assign damage state
    #                                       = 'Probablistic' to use lognormal fragilities to assign damage state
    #
    # OUTPUTS
    #    FDI         = np.array with all Floor Damage Index [0.0 <= FDI <= 1.0 all beams and columns in last damage state]
    #    DI          = Building Damage Indicator as the sum of the FDI [0 <= DI <= n_floors]
    #    DI_columns  = Building Damage Indicator accounting for beams only
    #    DI_beams    = Building Damage Indicator accounting for columns only
    #    columnBool  = 'True' if column damage took place (DI_columns/DI_beams > 0.00001)

    # Compute DI based on hinge demands
    n_stories, n_bays, _, _ = hinge_yield_rotation_positive.shape
    n_stories = n_stories - 1

    # Compute floor damage index for each floor
    FDI = np.zeros([n_stories + 1])
    FDIbeams = np.zeros([n_stories + 1])
    FDIcol = np.zeros([n_stories + 1])
    for floor_i in range(n_stories + 1):
        FDI[floor_i], FDIbeams[floor_i], FDIcol[floor_i] = compute_FDI(floor_i, peak_joint_pos, peak_joint_neg,
                                                                       hinge_yield_rotation_positive,
                                                                       hinge_cap_rotation_positive,
                                                                       hinge_yield_rotation_negative,
                                                                       hinge_cap_rotation_negative,
                                                                       inspector_type)

    # Computes the building's damage index as sum of all Floor Damage Indexes
    DI = sum(FDI)
    DI_columns = sum(FDIcol[1:]) # sum(FDIcol[1:]) ignores columns at first story since those yield very early
    DI_beams = sum(FDIbeams)

    if DI_columns / DI_beams > 0.00001:  # Not exactly 0 since there is a chance for at least one column damage
        columnBool = True
    else:
        columnBool = False

    return FDI, DI, DI_columns, DI_beams, columnBool


def design_spectra(T, Sms, Sm1):
    # Function to compute spectral acceleration Sa(T1) using ASCE7 code specification

    Tc = Sm1 / Sms

    Sa = np.zeros(len(T))
    for i in range(len(T)):
        Sa[i] = Sms
        if T[i] > Tc:
            Sa[i] = Sm1 / T[i]
    return Sa


def design_sa_avg(T1, Sms, Sm1):
    # Function to compute the Sa(T1)_avg using ASCE7 code specification and Eads et al. 2015 definition of Sa(T1)_avg

    T = np.linspace(0.2 * T1, 3 * T1, int(2.8 * T1 / 0.05))
    Sa = design_spectra(T, Sms, Sm1)

    Sa_avg = sp.stats.gmean(Sa)
    return Sa_avg


def sum_fdi_at_peak(idx_peak, fdi_i, n_floors_to_add):
    # Computes a damage indicators based on provided floor damage indexes summing n_floors around the story with
    # the peak story drift
    #
    # INPUTS
    #   idx_peak        = upper floor of the story with the peak story drift
    #   fdi_i           = np.array with the floor damage index per floor
    #   n_floors_to_add = the number of floors to add around the peak story drift
    #
    # OUTPUTS
    #   DI_sdr          = average damage index around peak story drift ratio
    #

    if idx_peak == 0:
        idx_peak += 1

    if n_floors_to_add == 1:
        a = fdi_i[idx_peak]
    else:
        if idx_peak <= np.ceil(n_floors_to_add / 2):
            a = fdi_i[0:n_floors_to_add]
        else:
            #             print(idx_peak - np.ceil(n_floors_to_add/2))
            #             print(idx_peak + n_floors_to_add - np.ceil(n_floors_to_add/2) - 1)
            #             print(fdi_i)
            a = fdi_i[int(idx_peak - np.ceil(n_floors_to_add / 2)):
                      int(idx_peak + n_floors_to_add - np.ceil(n_floors_to_add / 2) - 1)]

    return np.sum(a)/n_floors_to_add

def plot_hinges_prob(t, edp, joints_x, joints_y, plot_scale, peak_joint_pos, peak_joint_neg,
                     hinge_yield_rotation_positive,
                     hinge_cap_rotation_positive, hinge_yield_rotation_negative, hinge_cap_rotation_negative, ax):
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
    disp_t = np.insert(disp_t, 0, 0,
                       axis=0)  # add the hinge at column base# add zero displacement at base of the column
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

                # Define color for plot
                if yieldPos != 0:
                    # Compute expected di of current hinge
                    if hinge_loc == 1 or hinge_loc == 3:
                        isBeam = True
                        pos_fragilities = get_FEMAP58_fragility(isBeam, capPos)
                        neg_fragilities = get_FEMAP58_fragility(isBeam, capNeg)
                        di, _ = compute_di_prob(peakPos, peakNeg, pos_fragilities, neg_fragilities)
                    else:
                        isBeam = False
                        pos_fragilities = get_FEMAP58_fragility(isBeam, capPos)
                        neg_fragilities = get_FEMAP58_fragility(isBeam, capNeg)
                        di, _ = compute_di_prob(peakPos, peakNeg, pos_fragilities, neg_fragilities)

                    # Select color of current hinge based on di
                    if (di > 0.5) and (di <= 1.5):
                        # Damage state 1
                        joints_x_yield = np.append(joints_x_yield, curr_x)
                        joints_y_yield = np.append(joints_y_yield, curr_y)

                    elif (di > 1.5) and (di <= 2.5):
                        # Between capping/3 and capping*2/3
                        joints_x_capt = np.append(joints_x_capt, curr_x)
                        joints_y_capt = np.append(joints_y_capt, curr_y)

                    elif (di > 2.5):
                        # Between capping*2/3 and capping
                        #                         print('yieldPos='+str(yieldPos)+';capPos='+str(capPos))
                        #                         print('floor='+str(floor_i)+';column='+str(col_i)+';hinge='+str(hinge_loc))
                        joints_x_cap2t = np.append(joints_x_cap2t, curr_x)
                        joints_y_cap2t = np.append(joints_y_cap2t, curr_y)

    # Plot hinges with appropiate colors
    _ = ax.plot(joints_x_yield, joints_y_yield, 'o', color='m')
    _ = ax.plot(joints_x_capt, joints_y_capt, 'o', color='c')
    _ = ax.plot(joints_x_cap2t, joints_y_cap2t, 'o', color='r')


def get_max_edp_ratio(peak_joint_pos, peak_joint_neg, hinge_cap_rotation_positive, hinge_cap_rotation_negative):
    # Function to find the worst beam and column (the beam and column with the largest rotation demand)
    #
    # INPUTS
    #    peak_joint_pos                     = 4D array with the positive peak rotation in each hinge [floor_i,col_i,hinge_loc, 0]
    #                                         hinge_loc = 0 bottom; 1: right; 2: top; 3: left of joint at floor_i, col_i
    #    peak_joint_neg                     = 4D array with the negative peak rotation in each hinge [floor_i,col_i,hinge_loc, 0]
    #                                         hinge_loc = 0 bottom; 1: right; 2: top; 3: left of joint at floor_i, col_i
    #    hinge_cap_rotation_positive        = 4D array with the pos rotation at capping in each hinge [floor_i,col_i,hinge_loc, 0]
    #    hinge_cap_rotation_negative        = 4D array with the neg rotation at capping in each hinge [floor_i,col_i,hinge_loc, 0]
    #
    # OUTPUTS
    #    beam_max_rot_ratio                 = max_rot of the worst beam divided by its capping rotation
    #    col_max_rot_ratio                  = max_rot of the worst column divided by its capping rotation
    #

    # Basic geometrical data
    n_floors, n_columns, _, _ = hinge_cap_rotation_positive.shape

    beam_max_rot_ratio = 0
    col_max_rot_ratio = 0

    # Loop on all hinges
    for floor_i in range(n_floors):

        for col_i in range(n_columns):

            for hinge_loc in range(4):

                # Read rotation demand of current hinge
                peakPos = peak_joint_pos[floor_i, col_i, hinge_loc, 0]
                peakNeg = -peak_joint_neg[floor_i, col_i, hinge_loc, 0]

                # Read rotation capacity of current hinge
                capPos = hinge_cap_rotation_positive[floor_i, col_i, hinge_loc, 0]
                capNeg = -hinge_cap_rotation_negative[floor_i, col_i, hinge_loc, 0]

                # find worst element (di)
                if capPos != 0:
                    if hinge_loc == 1 or hinge_loc == 3:
                        isBeam = True
                        rot_ratio = max(peakPos / capPos, peakNeg / capNeg)

                        if rot_ratio > beam_max_rot_ratio:
                            beam_max_rot_ratio = rot_ratio
                    else:
                        isBeam = False
                        rot_ratio = max(peakPos / capPos, peakNeg / capNeg)

                        if rot_ratio > beam_max_rot_ratio:
                            col_max_rot_ratio = rot_ratio
    return beam_max_rot_ratio, col_max_rot_ratio


def compute_di_sim(peakPos, peakNeg, pos_fragilities, neg_fragilities):
    # Function to compute hinge damage index (di) accounting for fragilities with a monte-carlo approach
    # here we sample a unifirm random number per component, so we are assuming independende of damage states
    #
    # INPUTS
    #    peakPos          = float edp demand in positive direction
    #    peakNeg          = float edp demand in negative direction
    #    pos_fragilitieS  = DataFrame with damage fragilities for positive demand
    #    neg_fragilities  = DataFrame with damage fragilities for negative demand
    #
    # OUTPUTS
    #    di               = Simulated damage state
    #

    n_damage_states = len(pos_fragilities['median'])

    cumm_prob = np.zeros(n_damage_states)
    weights = np.zeros(n_damage_states + 1)

    # Generate uniform random number
    U = np.random.uniform()

    # Compute di for positive peak demand
    for ds_i in range(n_damage_states):
        median = pos_fragilities['median'][ds_i + 1]
        sigma_ln = pos_fragilities['beta'][ds_i + 1]
        cumm_prob[ds_i] = stats.lognorm.cdf(peakPos, sigma_ln, 0, median)

    for ds_i in range(n_damage_states + 1):
        if ds_i == 0:
            if U >= cumm_prob[ds_i]:
                di_pos = 0
        elif ds_i == n_damage_states:
            if U < cumm_prob[ds_i - 1]:
                di_pos = n_damage_states
        else:
            if U >= cumm_prob[ds_i] and U < cumm_prob[ds_i - 1]:
                di_pos = ds_i

    # Compute di for negative peak demand
    for ds_i in range(n_damage_states):
        median = neg_fragilities['median'][ds_i + 1]
        sigma_ln = neg_fragilities['beta'][ds_i + 1]
        cumm_prob[ds_i] = stats.lognorm.cdf(peakNeg, sigma_ln, 0, median)

    for ds_i in range(n_damage_states + 1):
        if ds_i == 0:
            if U >= cumm_prob[ds_i]:
                di_neg = 0
        elif ds_i == n_damage_states:
            if U < cumm_prob[ds_i - 1]:
                di_neg = n_damage_states
        else:
            if U >= cumm_prob[ds_i] and U < cumm_prob[ds_i - 1]:
                di_neg = ds_i

    return max(di_pos, di_neg)


def get_dsr_monte_carlo(peak_joint_pos, peak_joint_neg, hinge_cap_rotation_positive, hinge_cap_rotation_negative,
                        n_simulations):
    # get_dsr_monte_carlo computes the DAMAGE STATE RATIO of the beams and columns in a frame this is done without any
    # assumption of independence thus needs several simulations to account for possible correlatio.
    # This is by definition a probabilisitc approach
    #
    # INPUTS
    #    peak_joint_pos                     = 4D array with the positive peak rotation in each hinge [floor_i,col_i,hinge_loc, 0]
    #                                         hinge_loc = 0 bottom; 1: right; 2: top; 3: left of joint at floor_i, col_i
    #    peak_joint_neg                     = 4D array with the negative peak rotation in each hinge [floor_i,col_i,hinge_loc, 0]
    #                                         hinge_loc = 0 bottom; 1: right; 2: top; 3: left of joint at floor_i, col_i
    #    hinge_yield_rotation_positive      = 4D array with the pos yield rotation capacity in each hinge [floor_i,col_i,hinge_loc, 0]
    #    hinge_yield_rotation_negative      = 4D array with the neg yield rotation capacity in each hinge [floor_i,col_i,hinge_loc, 0]
    #    hinge_cap_rotation_positive        = 4D array with the pos rotation at capping in each hinge [floor_i,col_i,hinge_loc, 0]
    #    hinge_cap_rotation_negative        = 4D array with the neg rotation at capping in each hinge [floor_i,col_i,hinge_loc, 0]
    #    n_simulations                      = Number of simulations for montecarlo
    #
    # OUTPUTS
    #    dsr_beams     = np.array with the portion of beams at each damage state in the entire building
    #                    ex. dsr_beams = [0.7, 0.15, 0.05, 0.10]: 70% beams in DS0 (no damage), 15% in DS1, 5% in DS2, and 10% in DS3
    #    dsr_columns   = np.array with the portion of columns at each damage state in the entire building

    # Basic geometrical data
    n_floors, n_columns, _, _ = hinge_cap_rotation_positive.shape

    # Initialize array to store di of all hinges in the building
    beam_di = np.zeros([(n_columns - 1) * (n_floors - 1) * 2, n_simulations])
    column_di = np.zeros([n_columns * (n_floors - 1) * 2, n_simulations])

    # No. of fragility functions per component
    n_fragilities = 3

    # Inititialize damage index
    FDI = np.zeros(n_floors)
    DI = np.zeros(n_simulations)

    # Initialize damage state ratios per simulation
    dsr_beams = np.zeros([n_fragilities + 1, n_simulations])
    dsr_columns = np.zeros([n_fragilities + 1, n_simulations])

    for sim_i in range(n_simulations):

        # index for beam and columns arrays
        beam_hinge_i = 0
        column_hinge_i = 0

        # Loop on all hinges
        for floor_i in range(n_floors):

            for col_i in range(n_columns):

                for hinge_loc in range(4):
                    # Read rotation demand of current hinge
                    peakPos = peak_joint_pos[floor_i, col_i, hinge_loc, 0]
                    peakNeg = -peak_joint_neg[floor_i, col_i, hinge_loc, 0]

                    # Read rotation capacity of current hinge
                    capPos = hinge_cap_rotation_positive[floor_i, col_i, hinge_loc, 0]
                    capNeg = -hinge_cap_rotation_negative[floor_i, col_i, hinge_loc, 0]

                    # hinge damage index (di)
                    if capPos != 0:  # and peakPos != 0 and peakNeg != 0:
                        if hinge_loc == 1 or hinge_loc == 3:
                            pos_fragilities = get_FEMAP58_fragility(True, capPos)
                            neg_fragilities = get_FEMAP58_fragility(True, capNeg)
                            beam_di[beam_hinge_i, sim_i] = compute_di_sim(peakPos, peakNeg, pos_fragilities,
                                                                          neg_fragilities)
                            beam_hinge_i = beam_hinge_i + 1

                        else:
                            pos_fragilities = get_FEMAP58_fragility(False, capPos)
                            neg_fragilities = get_FEMAP58_fragility(False, capNeg)
                            column_di[column_hinge_i, sim_i] = compute_di_sim(peakPos, peakNeg, pos_fragilities,
                                                                              neg_fragilities)
                            column_hinge_i = column_hinge_i + 1

                            # Compute FDI
            beam_hinges_per_story = (n_columns - 1) * 2
            if floor_i == 0 or floor_i == n_floors:
                # at first and last floors only columns on one side
                column_hinges_per_story = n_columns
            else:
                # al all other floors columns on top and bottom
                column_hinges_per_story = n_columns * 2
            floor_beam_di = beam_di[(beam_hinge_i - beam_hinges_per_story):beam_hinge_i, sim_i]
            floor_column_di = column_di[(column_hinge_i - n_columns * 2):column_hinge_i, sim_i]

            FDI[floor_i] = (0.5 * sum(floor_beam_di) / n_fragilities / beam_hinges_per_story +
                            0.5 * sum(floor_column_di) / n_fragilities / column_hinges_per_story)

        # Get No. of beams and columns in each damage state
        for ds_i in range(n_fragilities):
            dsr_beams[ds_i, sim_i] = np.sum(beam_di[:, sim_i] == ds_i) / len(beam_di[:, sim_i])
            dsr_columns[ds_i, sim_i] = np.sum(column_di[:, sim_i] == ds_i) / len(column_di[:, sim_i])

        # Compute DI hinges
        DI[sim_i] = sum(FDI)

    # Mean damage state ratios for beams and columns
    mean_dsr_beams = np.zeros(n_fragilities + 1)
    mean_dsr_columns = np.zeros(n_fragilities + 1)

    for ds_i in range(n_fragilities):
        mean_dsr_beams[ds_i] = np.mean(dsr_beams[ds_i, :])
        mean_dsr_columns[ds_i] = np.mean(dsr_columns[ds_i, :])

    return mean_dsr_beams, mean_dsr_columns, np.mean(DI)


def get_dsr(peak_joint_pos, peak_joint_neg, hinge_yield_rotation_positive, hinge_yield_rotation_negative,
            hinge_cap_rotation_positive, hinge_cap_rotation_negative, inspector_type):
    # get_dsr computes the DAMAGE STATE RATIO of the beams and columns in a frame assuming that the damage states of every
    # component are mutually independent
    #
    # INPUTS
    #    peak_joint_pos                     = 4D array with the positive peak rotation in each hinge [floor_i,col_i,hinge_loc, 0]
    #                                         hinge_loc = 0 bottom; 1: right; 2: top; 3: left of joint at floor_i, col_i
    #    peak_joint_neg                     = 4D array with the negative peak rotation in each hinge [floor_i,col_i,hinge_loc, 0]
    #                                         hinge_loc = 0 bottom; 1: right; 2: top; 3: left of joint at floor_i, col_i
    #    hinge_yield_rotation_positive      = 4D array with the pos yield rotation capacity in each hinge [floor_i,col_i,hinge_loc, 0]
    #    hinge_yield_rotation_negative      = 4D array with the neg yield rotation capacity in each hinge [floor_i,col_i,hinge_loc, 0]
    #    hinge_cap_rotation_positive        = 4D array with the pos rotation at capping in each hinge [floor_i,col_i,hinge_loc, 0]
    #    hinge_cap_rotation_negative        = 4D array with the neg rotation at capping in each hinge [floor_i,col_i,hinge_loc, 0]
    #    inspector_type                     = 'Deterministic' to use deterministic fragilities to assign damage state
    #                                       = 'Probablistic' to use lognormal fragilities to assign damage state
    #
    # OUTPUTS
    #    dsr_beams     = np.array with the portion of beams at each damage state in the entire building
    #                    ex. dsr_beams = [0.7, 0.15, 0.05, 0.10]: 70% beams in DS0 (no damage), 15% in DS1, 5% in DS2, and 10% in DS3
    #    dsr_columns   = np.array with the portion of columns at each damage state in the entire building

    # Basic geometrical data
    n_floors, n_columns, _, _ = hinge_yield_rotation_positive.shape

    # Initialize array to store di of all hinges in the building
    beam_di = np.zeros([(n_columns - 1) * (n_floors - 1) * 2])
    column_di = np.zeros([n_columns * (n_floors - 1) * 2])
    if inspector_type == 'Probabilistic':
        weights_beams = np.zeros(
            [(n_columns - 1) * (n_floors - 1) * 2, 4])  # assuming 4 damage states (including no damage)
        weights_columns = np.zeros([n_columns * (n_floors - 1) * 2, 4])

    # index for beam and columns arrays
    beam_hinge_i = 0
    column_hinge_i = 0

    # Loop on all hinges
    for floor_i in range(n_floors):

        for col_i in range(n_columns):

            for hinge_loc in range(4):

                # Read rotation demand of current hinge
                peakPos = peak_joint_pos[floor_i, col_i, hinge_loc, 0]
                peakNeg = -peak_joint_neg[floor_i, col_i, hinge_loc, 0]

                # Read rotation capacity of current hinge
                yieldPos = hinge_yield_rotation_positive[floor_i, col_i, hinge_loc, 0]
                yieldNeg = -hinge_yield_rotation_negative[floor_i, col_i, hinge_loc, 0]
                capPos = hinge_cap_rotation_positive[floor_i, col_i, hinge_loc, 0]
                capNeg = -hinge_cap_rotation_negative[floor_i, col_i, hinge_loc, 0]

                # hinge damage index (di)
                if yieldPos != 0:  # and peakPos != 0 and peakNeg != 0:
                    if hinge_loc == 1 or hinge_loc == 3:
                        if inspector_type == 'Deterministic':
                            beam_di[beam_hinge_i] = compute_di_deterministic(True, peakPos, peakNeg, yieldPos, capPos,
                                                                             yieldNeg, capNeg)
                        else:
                            pos_fragilities = get_FEMAP58_fragility(True, capPos)
                            neg_fragilities = get_FEMAP58_fragility(True, capNeg)
                            beam_di[beam_hinge_i], weights_beams[beam_hinge_i, :] = compute_di_prob(peakPos, peakNeg, pos_fragilities, neg_fragilities)

                        beam_hinge_i = beam_hinge_i + 1

                    else:
                        if inspector_type == 'Deterministic':
                            column_di[column_hinge_i] = compute_di_deterministic(False, peakPos, peakNeg, yieldPos,
                                                                                 capPos, yieldNeg, capNeg)
                        else:
                            pos_fragilities = get_FEMAP58_fragility(False, capPos)
                            neg_fragilities = get_FEMAP58_fragility(False, capNeg)
                            column_di[column_hinge_i], weights_columns[column_hinge_i, :] = compute_di_prob(peakPos, peakNeg, pos_fragilities, neg_fragilities)

                        column_hinge_i = column_hinge_i + 1

    # No. of damage states per component
    if inspector_type == 'Deterministic':
        n_fragilities = 4
    else:
        n_fragilities = len(pos_fragilities)

    # Compute Damage State Ratios
    dsr_beams = np.zeros(n_fragilities + 1)
    dsr_columns = np.zeros(n_fragilities + 1)

    # Get No. of beams and columns in each damage state
    if inspector_type == 'Deterministic':
        for ds_i in range(n_fragilities + 1):
            dsr_beams[ds_i] = np.sum(beam_di == ds_i) / len(beam_di)
            dsr_columns[ds_i] = np.sum(column_di == ds_i) / len(column_di)
    else:
        for ds_i in range(n_fragilities + 1):
            dsr_beams[ds_i] = np.sum(weights_beams[:, ds_i]) / len(beam_di)
            dsr_columns[ds_i] = np.sum(weights_columns[:, ds_i]) / len(column_di)

    return dsr_beams, dsr_columns

def fun_to_fit3Lin(x, k_max, a_min, a, b1, b2):
    y = np.zeros(x.size)

    for i in range(x.size):
        if x[i] < a_min:
            y[i] = k_max
        elif x[i] < a:
            y[i] = k_max + b1 * (x[i] - a_min)
        else:
            y[i] = k_max + b1 * (a - a_min) + b2 * (x[i] - a)

    return y


def least_squared_deviations_3Lin(param, x, y):
    # loss function for the parameters in param list

    cost = np.sum((fun_to_fit3Lin(x, param[0], param[1], param[2], param[3], param[4]) - y) ** 2)

    #     print(param)
    #     print(cost)
    return cost


def grad_vector(param, x, y):
    # Gradient vector of the loss function for the parameters in param list

    grad = np.zeros(param.size)

    y_pred = fun_to_fit3Lin(x, param[0], param[1], param[2], param[3], param[4])

    y_seg1 = y[np.logical_and(x > param[1], x <= param[2])]
    y_seg2 = y[x > param[2]]

    y_pred_seg1 = y_pred[np.logical_and(x > param[1], x <= param[2])]
    y_pred_seg2 = y_pred[x > param[2]]

    x_seg1 = x[np.logical_and(x > param[1], x <= param[2])]
    x_seg2 = x[x > param[2]]

    grad[0] = 2 * np.sum(y_pred - y)  # derivative with respect to k_max
    grad[1] = 2 * (np.sum((y_pred_seg1 - y_seg1) * -param[3]) + np.sum(
        (y_pred_seg2 - y_seg2) * -param[3]))  # derivative with respect to a_min
    grad[2] = 2 * (np.sum((y_pred_seg2 - y_seg2) * -param[4]))  # derivative with respect to a
    grad[3] = 2 * (np.matmul(np.transpose(y_pred_seg1 - y_seg1), x_seg1 - param[1]) +
                   np.sum((y_pred_seg2 - y_seg2) * (param[2] - param[1])))  # derivative with respect to b1
    grad[4] = 2 * np.matmul(np.transpose(y_pred_seg2 - y_seg2), x_seg2 - param[2])  # derivative with respect to b2

    return grad


def fitPieceWiseFunc3LinLS_oneGuess(DIname, x, y, x_min, x_data_space):
    # Fit the piecewise linear function minimizing the sum of squared deviation. Starts by using a tri-linear function
    # with first segment horizontal in 1.0. If the x-value for the first kink is very close to the minimum value of the
    # vector x, there is a possible cluster of data at x = 0 and changes the function for a bi-linear function that does
    # not need to start in 1 but at a lower value.
    #
    # INPUTS
    #    x = 1D np.array with independent variable (damage indicator)
    #    y = 1D np/array with dependent variable (reduction in collapse capacity)
    #    x_min = Minimum value to consider for damage indicators (very important if using log space, the function replaces
    #            zeros by this small number)
    #    x_data_space = Space for fitting 'Linear' no change in original data
    #                    'Log' takes the natural log
    #                    'log_std' takes the natural log and then standarized the data
    # OUTPUTS
    #    parameters_DI = pd.DataFrame witht the information of the fit
    #                   'Minimum limit'
    #                   'Threshold limit'
    #                   'Slope 1'
    #                   'Slope 2'
    #                   'k_max'
    #                   'std_residuals'
    #                   'Residuals'
    #                   'x_min'
    #                   'di for plot'
    #                   'k for plot'

    ##### Transform data if requiered #####
    x = x.flatten()

    x[x <= x_min] = x_min

    if x_data_space == "log":
        x = np.log(x)
    elif x_data_space == "std":
        mu_x = np.mean(x)
        std_x = np.std(x)
        x = (x - mu_x) / std_x
    elif x_data_space == "log_std":
        x = np.log(x)
        mu_lnx = np.mean(x)
        std_lnx = np.std(x)
        x = (x - np.mean(x)) / np.std(x)

    ##### Initial guess #####
    a0 = np.mean(x)
    a0_min = np.mean([min(x), a0])

    x1 = x[x < a0] - min(x)
    y1 = y[x < a0] - 1
    b1_0 = np.sum(x1.dot(y1)) / np.sum(x1.dot(x1))

    x2 = x[x >= a0] - a0
    y2 = y[x >= a0] - (1 + b1_0 * a0)
    b2_0 = np.sum(x2.dot(y2)) / np.sum(x2.dot(x2))

    k_max_0 = 1.0

    param0 = [k_max_0, a0_min, a0, b1_0, b2_0]

    #     print(param0)

    ##### Optimization #####

    # Constraint that a_min > a_
    def cons(param):
        return (param[2] - param[1])

    # Boundaries for parameters
    bnds = ((0.9, 1.1), (min(x), max(x)), (min(x), max(x)), (None, None), (None, None))

    # Perform optimization
    output = optimize.minimize(least_squared_deviations_3Lin, param0, method='SLSQP', jac=grad_vector, args=(x, y),
                               bounds=bnds, constraints={"fun": cons, "type": "ineq"})
    k_max = output.x[0]
    amin_opt = output.x[1]
    a_opt = output.x[2]
    b1_opt = output.x[3]
    b2_opt = output.x[4]
    print(output.fun)

    ##### Return to original space vector for plots #####
    a_plot = np.array([min(x), amin_opt, a_opt, max(x)])
    k_plot = fun_to_fit3Lin(a_plot, k_max, amin_opt, a_opt, b1_opt, b2_opt)

    if x_data_space == "log":
        a_plot = np.exp(a_plot)
    elif x_data_space == "std":
        a_plot = a_plot * std_x + mu_x
    elif x_data_space == "log_std":
        a_plot = np.exp(a_plot * std_lnx + mu_lnx)

    ##### Efficiency #####
    ypred = fun_to_fit3Lin(x, k_max, amin_opt, a_opt, b1_opt, b2_opt)

    residual = y - ypred
    std_res = np.std(residual)

    columns = ['Minimum limit', 'Threshold limit', 'Slope 1', 'Slope 2', 'k_max', 'std_residuals', 'Residuals', 'x_min',
               'di for plot', 'k for plot']
    parameters_DI = pd.DataFrame([[amin_opt, a_opt, b1_opt, b2_opt, k_max, std_res, residual, x_min, a_plot, k_plot]],
                                 columns=columns, index=[DIname])

    return parameters_DI



def fitPieceWiseFunc3LinLS(DIname, x, y, x_min, x_data_space, N_seeds):
    # Fit the piecewise linear function minimizing the sum of squared deviation. Starts by using a tri-linear function
    # with first segment horizontal in 1.0. If the x-value for the first kink is very close to the minimum value of the
    # vector x, there is a possible cluster of data at x = 0 and changes the function for a bi-linear function that does
    # not need to start in 1 but at a lower value.
    #
    # INPUTS
    #    x = 1D np.array with independent variable (damage indicator)
    #    y = 1D np/array with dependent variable (reduction in collapse capacity)
    #    x_min = Minimum value to consider for damage indicators (very important if using log space, the function replaces
    #            zeros by this small number)
    #    x_data_space = Space for fitting 'Linear' no change in original data
    #                    'Log' takes the natural log
    #                    'log_std' takes the natural log and then standarized the data
    #    N_seeds = Number of random initial guesses of parameters
    #
    # OUTPUTS
    #    parameters_DI = pd.DataFrame witht the information of the fit
    #                   'Minimum limit'
    #                   'Threshold limit'
    #                   'Slope 1'
    #                   'Slope 2'
    #                   'k_max'
    #                   'std_residuals'
    #                   'Residuals'
    #                   'x_min'
    #                   'di for plot'
    #                   'k for plot'

    ##### Transform data if requiered #####
    x = x.flatten()

    x[x <= x_min] = x_min

    if x_data_space == "log":
        x = np.log(x)
    elif x_data_space == "std":
        mu_x = np.mean(x)
        std_x = np.std(x)
        x = (x - mu_x) / std_x
    elif x_data_space == "log_std":
        x = np.log(x)
        mu_lnx = np.mean(x)
        std_lnx = np.std(x)
        x = (x - np.mean(x)) / np.std(x)

    k_max_trial = np.zeros(N_seeds)
    amin_trial = np.zeros(N_seeds)
    a_trial = np.zeros(N_seeds)
    b1_trial = np.zeros(N_seeds)
    b2_trial = np.zeros(N_seeds)
    MSE_trial = np.empty(N_seeds)
    MSE_trial[:] = np.NaN

    for i_trial in range(N_seeds):
        while np.isnan(MSE_trial[i_trial]):
            ##### Initial guess #####
            a0 = random.uniform(min(x), max(x))
            a0_min = random.uniform(min(x), a0)

            k_max_0 = 1.0

            x1 = x[x < a0] - min(x)
            y1 = y[x < a0] - 1
            b1_0 = np.sum(x1.dot(y1)) / np.sum(x1.dot(x1))

            x2 = x[x >= a0] - a0
            y2 = y[x >= a0] - (1 + b1_0 * a0)
            b2_0 = np.sum(x2.dot(y2)) / np.sum(x2.dot(x2))

            param0 = [k_max_0, a0_min, a0, b1_0, b2_0]

            ##### Optimization #####
            # Constraint that a_min > a_
            def cons(param):
                return (param[2] - param[1])

            # Boundaries for parameters
            bnds = ((0.95, 1.05), (min(x), max(x)), (min(x), max(x)), (None, None), (None, None))

            # Perform optimization
            output = optimize.minimize(least_squared_deviations_3Lin, param0, method='SLSQP', jac=grad_vector, args=(x, y),
                                       bounds=bnds, constraints={"fun": cons, "type": "ineq"})
            k_max_trial[i_trial] = output.x[0]
            amin_trial[i_trial] = output.x[1]
            a_trial[i_trial] = output.x[2]
            b1_trial[i_trial] = output.x[3]
            b2_trial[i_trial] = output.x[4]
            MSE_trial[i_trial] = output.fun

    #### Select best solution ####
    i_opt = np.argwhere(MSE_trial == min(MSE_trial))
    # print(MSE_trial)
    # print(MSE_trial[i_opt])

    k_max = float(k_max_trial[i_opt[0]])
    amin_opt = float(amin_trial[i_opt[0]])
    a_opt = float(a_trial[i_opt[0]])
    b1_opt = float(b1_trial[i_opt[0]])
    b2_opt = float(b2_trial[i_opt[0]])

    ##### Return to original space vector for plots #####
    a_plot = np.array([min(x), amin_opt, a_opt, max(x)])
    k_plot = fun_to_fit3Lin(a_plot, k_max, amin_opt, a_opt, b1_opt, b2_opt)

    if x_data_space == "log":
        a_plot = np.exp(a_plot)
    elif x_data_space == "std":
        a_plot = a_plot * std_x + mu_x
    elif x_data_space == "log_std":
        a_plot = np.exp(a_plot * std_lnx + mu_lnx)

    ##### Efficiency #####
    ypred = fun_to_fit3Lin(x, k_max, amin_opt, a_opt, b1_opt, b2_opt)

    residual = y - ypred
    std_res = np.std(residual)

    columns = ['Minimum limit', 'Threshold limit', 'Slope 1', 'Slope 2', 'k_max', 'std_residuals', 'Residuals', 'x_min',
               'di for plot', 'k for plot']
    parameters_DI = pd.DataFrame([[amin_opt, a_opt, b1_opt, b2_opt, k_max, std_res, residual, x_min, a_plot, k_plot]],
                                 columns=columns, index=[DIname])

    return parameters_DI


def predictPieceWiseFunc3LinLS(x, parameters, x_data_space):
    # Uses the fitted tri-linear (or bilienar) function to estimate the values of y
    #
    # INPUTS
    #    x = 1D np.array with independent variable (damage indicator)
    #    parameters   = pd.DataFrame with the information of the fitted tri-linear function in the same x_data_space
    #    x_data_space = Space for fitting 'Linear' no change in original data
    #                    'Log' takes the natural log
    #                    'log_std' takes the natural log and then standarized the data
    #
    # OUTPUTS
    #    y = prediction

    x = x.flatten()

    a_opt = parameters["Threshold limit"][0]
    amin_opt = parameters["Minimum limit"][0]
    b1_opt = parameters["Slope 1"][0]
    b2_opt = parameters["Slope 2"][0]
    x_min = parameters["x_min"][0]
    k_max = parameters["k_max"][0]

    x[x <= x_min] = x_min

    if x_data_space == "log":
        x = np.log(x)
    elif x_data_space == "std":
        mu_x = np.mean(x)
        std_x = np.std(x)
        x = (x - mu_x) / std_x
    elif x_data_space == "log_std":
        x = np.log(x)
        mu_lnx = np.mean(x)
        std_lnx = np.std(x)
        x = (x - np.mean(x)) / np.std(x)

    # Initial guess
    y = fun_to_fit3Lin(x, k_max, amin_opt, a_opt, b1_opt, b2_opt)

    return y

##### 2-piece wise linear function #####
def fun_to_fit2Lin(x, k_max, a_min, a, b1, b2):
    y = np.zeros(x.size)

    for i in range(x.size):
        if x[i] < a_min:
            y[i] = k_max
        elif x[i] < a:
            y[i] = k_max + b1 * (x[i] - a_min)
        else:
            y[i] = k_max + b1 * (a - a_min) + b2 * (x[i] - a)

    return y


def least_squared_deviations_2Lin(param, x, y, a_min):
    # loss function for the parameters in param list

    cost = np.sum((fun_to_fit2Lin(x, param[0], a_min, param[1], param[2], param[3]) - y) ** 2)

    return cost


def grad_vector_2Lin(param, x, y, a_min):
    # Gradient vector of the loss function for the parameters in param list

    grad = np.zeros(param.size)

    y_pred = fun_to_fit2Lin(x, param[0], a_min, param[1], param[2], param[3])

    y_seg1 = y[np.logical_and(x > a_min, x <= param[1])]
    y_seg2 = y[x > param[1]]

    y_pred_seg1 = y_pred[np.logical_and(x > a_min, x <= param[1])]
    y_pred_seg2 = y_pred[x > param[1]]

    x_seg1 = x[np.logical_and(x > a_min, x <= param[1])]
    x_seg2 = x[x > param[1]]

    grad[0] = 2 * np.sum(y_pred - y)  # derivative with respect to k_max
    grad[1] = 2 * (np.sum((y_pred_seg2 - y_seg2) * -param[3]))  # derivative with respect to a
    grad[2] = 2 * (np.matmul(np.transpose(y_pred_seg1 - y_seg1), x_seg1 - a_min) +
                   np.sum((y_pred_seg2 - y_seg2) * (param[1] - a_min)))  # derivative with respect to b1
    grad[3] = 2 * np.matmul(np.transpose(y_pred_seg2 - y_seg2), x_seg2 - param[1])  # derivative with respect to b2

    return grad


def fitPieceWiseFunc2LinLS_oneGuess(DIname, x, y, a_min, x_data_space):
    # Fit the piecewise linear function minimizing the sum of squared deviation. Starts by using a tri-linear function
    # with first segment horizontal in 1.0. If the x-value for the first kink is very close to the minimum value of the
    # vector x, there is a possible cluster of data at x = 0 and changes the function for a bi-linear function that does
    # not need to start in 1 but at a lower value.
    #
    # INPUTS
    #    x = 1D np.array with independent variable (damage indicator)
    #    y = 1D np/array with dependent variable (reduction in collapse capacity)
    #    x_min = Minimum value to consider for damage indicators (very important if using log space, the function replaces
    #            zeros by this small number)
    #    x_data_space = Space for fitting 'Linear' no change in original data
    #                    'Log' takes the natural log
    #                    'log_std' takes the natural log and then standarized the data
    # OUTPUTS
    #    parameters_DI = pd.DataFrame witht the information of the fit
    #                   'Minimum limit'
    #                   'Threshold limit'
    #                   'Slope 1'
    #                   'Slope 2'
    #                   'k_max'
    #                   'std_residuals'
    #                   'Residuals'
    #                   'x_min'
    #                   'di for plot'
    #                   'k for plot'

    ##### Transform data if requiered #####
    x = x.flatten()

    x[x <= 1e-8] = 1e-8

    if x_data_space == "log":
        a_min = np.log(a_min)
        x = np.log(x)
    elif x_data_space == "std":
        mu_x = np.mean(x)
        std_x = np.std(x)
        x = (x - mu_x) / std_x
        a_min = (np.log(a_min) - mu_x) / std_x
    elif x_data_space == "log_std":
        x = np.log(x)
        mu_lnx = np.mean(x)
        std_lnx = np.std(x)
        x = (x - mu_lnx) / std_lnx
        a_min = (a_min - np.mean(x)) / np.std(x)

    ##### Initial guess #####
    a0 = np.mean(x[x > a_min])

    x1 = x[x < a0] - a_min
    y1 = y[x < a0] - 1
    b1_0 = np.sum(x1.dot(y1)) / np.sum(x1.dot(x1))

    x2 = x[x >= a0] - a0
    y2 = y[x >= a0] - (1 + b1_0 * a0)
    b2_0 = np.sum(x2.dot(y2)) / np.sum(x2.dot(x2))

    k_max_0 = 1.0

    param0 = [k_max_0, a0, b1_0, b2_0]

    ##### Optimization #####

    # Constraint that a_min > a_
    def cons(param):
        return (param[1] - a_min)

    # Boundaries for parameters

    bnds = ((0.9, 1.1), (min(x), max(x)), (None, None), (None, None))

    # Perform optimization
    output = optimize.minimize(least_squared_deviations_2Lin, param0, method='SLSQP', jac=grad_vector_2Lin, args=(x, y, a_min),
                               bounds=bnds, constraints={"fun": cons, "type": "ineq"})
    k_max = output.x[0]
    a_opt = output.x[1]
    b1_opt = output.x[2]
    b2_opt = output.x[3]

    ##### Return to original space vector for plots #####
    a_plot = np.array([min(x), a_min, a_opt, max(x)])
    k_plot = fun_to_fit2Lin(a_plot, k_max, a_min, a_opt, b1_opt, b2_opt)

    if x_data_space == "log":
        a_plot = np.exp(a_plot)
    elif x_data_space == "std":
        a_plot = a_plot * std_x + mu_x
    elif x_data_space == "log_std":
        a_plot = np.exp(a_plot * std_lnx + mu_lnx)

    ##### Efficiency #####
    ypred = fun_to_fit2Lin(x, k_max, a_min, a_opt, b1_opt, b2_opt)

    residual = y - ypred
    std_res = np.std(residual)

    columns = ['Minimum limit', 'Threshold limit', 'Slope 1', 'Slope 2', 'k_max', 'std_residuals', 'Residuals', 'x_min',
               'di for plot', 'k for plot']
    parameters_DI = pd.DataFrame([[a_min, a_opt, b1_opt, b2_opt, k_max, std_res, residual, a_min, a_plot, k_plot]],
                                 columns=columns, index=[DIname])

    return parameters_DI



def plotDIvsk3Lin(x, y, x_limits, parameters, y_label, x_data_space, plot_i, color_specs, color_line, building_title):
    # Plot several damage indicators vs k including fitted function to observe their behavior
    #
    # INPUTS
    #     x            = 1D np.array with the values of the damage indicator to plot
    #     y            = 1D np.array with the values of the reduction in collapse capacity to plot
    #     x_limits     = [x_min, x_max] to plot
    #                    if 999 takes extremes of vector x for plot limits
    #                    if vector with one entry only, use it as minimum limit
    #     parameters   = pd.DataFrame with the information of the fitted tri-linear function in the same x_data_space
    #     y_label      = label for vertical axis
    #     x_data_space = Space for plotting 'Linear'
    #                    'Log'
    #                    'log_std'
    #     plot_i       = order for subplot (1 to 6)
    #     color_specs  = 1D np.array with color RBG for markers in scatter plot
    #     color_line   = 1D np.array with color RBG for regression line

    a_opt = parameters["Threshold limit"][0]
    amin_opt = parameters["Minimum limit"][0]
    b1_opt = parameters["Slope 1"][0]
    b2_opt = parameters["Slope 2"][0]
    a_plot = parameters["di for plot"][0]
    k_plot = parameters["k for plot"][0]
    x_min = parameters["x_min"][0]

    # Only plot data in the window requested of x
    if x_limits == 999:
        x_limits = [x_min, max(x)]
    elif len(x_limits) == 1:
        x_limits = [x_limits[0], max(x)]

    x[x <= x_min] = x_min

    x_label = parameters.index.values[0]

    if plot_i >= 0:
        # Locate the subplot
        if plot_i <= 2:
            ax = plt.subplot2grid((3, 3), (0, plot_i), rowspan=1, colspan=1)
        elif plot_i <= 5:
            ax = plt.subplot2grid((3, 3), (1, plot_i - 3), rowspan=1, colspan=1)
        else:
            ax = plt.subplot2grid((3, 3), (2, plot_i - 6), rowspan=1, colspan=1)

    # Plot analysis results
    ax = plt.scatter(x, y, s=10, color=color_specs, alpha=0.3)

    # Add reference line at Dmin and a
    ax = plt.plot([a_plot[2], a_plot[2]], [0, 1], linestyle='--', Color='gray', label='_nolegend_')
    ax = plt.plot([a_plot[1], a_plot[1]], [0, 1], linestyle='--', Color='gray', label='_nolegend_')

    # Add reference line in no reduction of collapse capacity
    ax = plt.plot(x_limits, [1, 1], linestyle='--', Color='gray', label='_nolegend_')

    # Add tri-piecewise linear regression of the data
    ax = plt.plot(a_plot, k_plot, linewidth=2.0, color=color_line, label='_nolegend_')

    # Format the plot and add title with the threshold
    ax = plt.xlabel(x_label)
    ax = plt.ylabel(y_label)
    ax = plt.xlim(x_limits)
    ax = plt.ylim([0.2, 1.1])
    ax = plt.grid(which='major', alpha=0.3)

    if x_data_space == "log_std" or x_data_space == "log":
        ax = plt.xscale('log')
        ax = plt.gca()
        _ = ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda y, _: '{:.2g}'.format(y)))

    # Title
    if parameters.index.values[0] == "DI hinges" or parameters.index.values[0] == "$FDI_{peak}$" or \
            parameters.index.values[0][0:3] == "$DI" or parameters.index.values[0][0:3] == "$Sa":
        title_text = "$a_{1}$ = " + str(round(a_plot[1], 4)) + "\n $a_{2}$ = " + str(round(a_plot[2], 2))
    else:
        title_text = "$a_{1}$ = " + str(round(a_plot[1] * 100, 2)) + "% \n $a_{2}$ = " + str(round(a_plot[2] * 100, 2)) + "%"
    if building_title != 0:
        title_text = building_title + "\n" + title_text

    ax = plt.title(title_text)

def plotResiduals(DI_name, x, parameters, y_label, x_data_space, plot_i):
    # Plot the residuals of
    #
    # INPUTS
    #     DI_name      = string with the name of the damage indicator to plot (same name as row in parameters DataFrame)
    #     x            = 1D np.array with the values of the damage indicator to plot
    #     parameters   = pd.DataFrame with the information of the fitted tri-linear function in the same x_data_space
    #     y_label      = label for vertical axis
    #     x_data_space = Space for plotting 'Linear'
    #                    'Log'
    #                    'log_std'
    #     plot_i       = order for subplot (1 to 6)

    # Retrieve residuals
    residuals = parameters["Residuals"][0]
    x_min = parameters["x_min"][0]

    # Only x greater than the minimum used in the fit
    x[x <= x_min] = x_min
    x_limits = [x_min, max(x)]

    x_label = parameters.index.values[0]

    # Choose appropiate subplot
    if plot_i <= 2:
        ax = plt.subplot2grid((3, 3), (0, plot_i), rowspan=1, colspan=1)
    elif plot_i <= 5:
        ax = plt.subplot2grid((3, 3), (1, plot_i - 3), rowspan=1, colspan=1)
    else:
        ax = plt.subplot2grid((3, 3), (2, plot_i - 6), rowspan=1, colspan=1)

    # Plot residuals
    _ = plt.scatter(x, residuals, s=10)

    # Format figure
    _ = plt.plot(x_limits, [0, 0], linestyle='--', Color='gray')
    _ = plt.ylabel(y_label)
    _ = plt.xlabel(x_label)


# Efficiency (Cross validation)
def cross_validation(x, y, x_min, x_transformation, N_seeds, splits=2, repeats=10):
    # cross_validation computes the mean squared error (MSE) of the fitted model by using k-fold cross validation
    # dividing the data (x) in splits groups to fit with one groups and compute the error with the other group
    #
    # INPUTS
    #     x = 1D np.array with independent variable (damage indicator)
    #     y = 1D np/array with dependent variable (reduction in collapse capacity)
    #     x_min = Minimum value to consider for damage indicators (very important if using log space, the function replaces
    #             zeros by this small number)
    #     x_transformation =  Space for fitting 'Linear' no change in original data
    #                    'Log' takes the natural log
    #                    'log_std' takes the natural log and then standarized the data
    #     splits           = number of segments to divide the dataset
    #     repeats          = number of repetitions of fitting-testing
    #
    # OUTPUT
    #     MAE = 1D np.array with the mean absolute error for all the repetitions

    MAE = np.zeros(repeats * splits)
    rkf = RepeatedKFold(n_splits=splits, n_repeats=repeats)
    j = 0
    for train_index, test_index in rkf.split(x):
        x_train, x_test = x[train_index], x[test_index]
        y_train, y_test = y[train_index], y[test_index]

        parameters = fitPieceWiseFunc3LinLS("trial", x_train, y_train, x_min, x_transformation, N_seeds)
#         parameters = fitPieceWiseFunc3LinLS_oneGuess("trial", x_train, y_train, x_min, x_transformation)

        y_pred = predictPieceWiseFunc3LinLS(x_test, parameters, x_transformation)

        MAE[j] = 1 / len(x_test) * np.sum(np.abs(y_pred - y_test))
        j = j + 1

    print('Done')
    return MAE

def sufficiencyPlot(gm_feature, parameters, x_label, y_label, plot_i):
    # Computes the p-value of the null-hypothesis that the coefficient of the linear regression is zero and plots the
    # residuals of the reduction in collapse capacity based on a property of the mainshock that caused the damage
    #
    # INPUTS
    #    gm_feature   = 1D np.array with the property of the mainshock to check sufficiency against
    #    parameters   = pd.DataFrame with the information of the fitted tri-linear function in the same x_data_space
    #    x_label      = label for horizontal axis
    #    y_label      = label for vertical axis
    #    plot_i       = order for subplot (0 to 5) use negative number to avoid plot

    # Retrieve residuals and x label
    residuals = parameters["Residuals"][0]
    DI_name = parameters.index.values[0]

    # Linear regression
    slope, intercept, _, p_value, _ = stats.linregress(gm_feature, residuals)

    if plot_i >= 0:
        # Choose appropiate subplot
        if plot_i <= 2:
            ax = plt.subplot2grid((3, 3), (0, plot_i), rowspan=1, colspan=1)
        elif plot_i <= 5:
            ax = plt.subplot2grid((3, 3), (1, plot_i - 3), rowspan=1, colspan=1)
        else:
            ax = plt.subplot2grid((3, 3), (2, plot_i - 6), rowspan=1, colspan=1)

        # Plot residuals vs gm feature to evaluate sufficiency
        ax = plt.scatter(gm_feature, residuals, alpha=0.3)

        # Add linear regression line
        x = np.linspace(0, max(gm_feature))
        resPred = intercept + slope * x

        # Figure format
        ax = plt.plot([0, max(gm_feature)], [0, 0], linestyle='--', Color='gray')
        ax = plt.ylabel(y_label)
        ax = plt.xlabel(x_label)
        ax = plt.plot(x, resPred, color="red")
        ax = plt.title(DI_name + "\n p-value = " + str(np.round(p_value, 4)))

    return p_value


def evaluate_thresholds(di_condition_1, threshold_1, di_condition_2, threshold_2, k, k_limit=0.8):
    # Counts the amount of false red tags, false green tag, and average fraction of correctly tagged buildings
    # using 2 damage indicators with their threshold and setting a limit in kappa for true safe vs un-safe
    #
    # INPUT
    #   di_condition_1  = np.array with the values of damage indicator 1 for each building damage instance
    #   threshold_1     = threshold to tag building damage instance per damage indicator 1
    #   di_condition_2  = np.array with the values of damage indicator 2 for each building damage instance
    #   threshold_2     = threshold to tag building damage instance per damage indicator 2
    #   k               = np.array of the reduction in collapse capacity for each building damage instance
    #   k_limit         = 0.8 (Limit to define safe vs un-safe buildings)
    #
    # OUTPUT
    #   fraction_FP     = false green tags / true green tags
    #   fraction_FN     = false red tags / true red tags
    #   accuracy        = average of the fraction of correctly tagged red and green building damage instances
    #

    # Compute safe and unsafe cases based on kappa
    unsafe_cases = np.count_nonzero(k < k_limit)
    safe_cases = np.count_nonzero(k >= k_limit)

    # Factors to remove dependency of number of examples in each class
    if max(safe_cases, unsafe_cases) == unsafe_cases:
        unsafe_factor = 1
        safe_factor = unsafe_cases / safe_cases
    else:
        unsafe_factor = safe_cases / unsafe_cases
        safe_factor = 1

    # Compute safe and unsafe cases predictions based on damage indicator conditions
    safe_pediction = np.count_nonzero([(di_condition_1 <= threshold_1) & (di_condition_2 <= threshold_2)])
    unsafe_pediction = np.count_nonzero([(di_condition_1 > threshold_1) | (di_condition_2 > threshold_2)])

    # Statistical measures
    FP = np.count_nonzero([(di_condition_1 <= threshold_1) & (di_condition_2 <= threshold_2) & (k < k_limit)])
    FN = np.count_nonzero([((di_condition_1 > threshold_1) | (di_condition_2 > threshold_2)) & (k > k_limit)])

    if safe_pediction != 0:
        fraction_FP = FP / unsafe_cases
    else:
        fraction_FP = 0

    if unsafe_pediction != 0:
        fraction_FN = FN / safe_cases
    else:
        fraction_FN = 0

    accuracy = ((unsafe_cases - FP) * unsafe_factor + (safe_cases - FN) * safe_factor) / (
                safe_cases * safe_factor + unsafe_cases * unsafe_factor)

    return fraction_FP, fraction_FN, accuracy


def evaluate_one_threshold(di_condition_1, threshold_1, k, k_limit=0.8):
    # Counts the amount of false red tags, false green tag, and average fraction of correctly tagged buildings
    # using 1 damage indicator with their threshold and setting a limit in kappa for true safe vs un-safe
    #
    # INPUT
    #   di_condition_1  = np.array with the values of damage indicator 1 for each building damage instance
    #   threshold_1     = threshold to tag building damage instance per damage indicator 1
    #   k               = np.array of the reduction in collapse capacity for each building damage instance
    #   k_limit         = 0.8 (Limit to define safe vs un-safe buildings)
    #
    # OUTPUT
    #   fraction_FP     = false green tags / true green tags
    #   fraction_FN     = false red tags / true red tags
    #   accuracy        = average of the fraction of correctly tagged red and green building damage instances
    #

    # Compute safe and unsafe cases based on kappa
    unsafe_cases = np.count_nonzero(k < k_limit)
    safe_cases = np.count_nonzero(k >= k_limit)

    # Factors to remove dependency of number of examples in each class
    if max(safe_cases, unsafe_cases) == unsafe_cases:
        unsafe_factor = 1
        safe_factor = unsafe_cases / safe_cases
    else:
        unsafe_factor = safe_cases / unsafe_cases
        safe_factor = 1

    # Compute safe and unsafe cases predictions based on damage indicator conditions
    safe_pediction = np.count_nonzero([(di_condition_1 <= threshold_1)])
    unsafe_pediction = np.count_nonzero([(di_condition_1 > threshold_1)])

    # Statistical measures
    FP = np.count_nonzero([(di_condition_1 <= threshold_1) & (k < k_limit)])
    FN = np.count_nonzero([(di_condition_1 > threshold_1) & (k > k_limit)])

    if safe_pediction != 0:
        fraction_FP = FP / unsafe_cases # Fraction of wrong greens per correct red tag
        #fraction_FP = FP / safe_prediction # Fraction of wrong greens per predicted green tags
    else:
        fraction_FP = 0

    if unsafe_pediction != 0:
        fraction_FN = FN / safe_cases
        #fraction_FN = FN / unsafe_prediction
    else:
        fraction_FN = 0

    accuracy = ((unsafe_cases - FP) * unsafe_factor + (safe_cases - FN) * safe_factor) / (
            safe_cases * safe_factor + unsafe_cases * unsafe_factor)

    return fraction_FP, fraction_FN, accuracy