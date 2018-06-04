import json
import numpy as np

N_TEST_TRIALS_PER_BLOCK = 32

# expected distance from random selection strategy
RAND_DISTANCE_1D = .25
RAND_DISTANCE_2D = .235

MED_SAMPLE_DIST_1D = .25
MED_SAMPLE_DIST_2D = .217

DATADIR = '../../data_raw/experiment1'

DATA = {}
COND_LABELS = ['1D-DIAL', '1D-RECT', '2D-DIAL', '2D-RECT']

SEL_CRITERION = .05


def data(sid):
    if sid in DATA:
        return DATA[sid]
    else:
        with open('%s/%s.json' % (DATADIR, sid), 'r') as f:
            DATA[sid] = json.load(f)
        return DATA[sid]

def condition(sid):
    return filter(lambda d: d[0]=='condition', data(sid))[0][1]

def rule_condition(sid):
    return filter(lambda d: d[0]=='rule_cond', data(sid))[0][1]

def stim_condition(sid):
    return filter(lambda d: d[0]=='stim_cond', data(sid))[0][1]

def rule_counter(sid):
    return filter(lambda d: d[0]=='counterbalance', data(sid))[0][1]

def dim_mapping(sid):
    return filter(lambda d: d[0]=='dim_mapping', data(sid))[0][1]

def rule_type(sid):
    cond = rule_condition(sid)
    counter = rule_counter(sid)
    m = dim_mapping(sid)
    stim = stim_condition(sid)

    if cond=='rb':
        # 1D
        if (counter, m) in [(0,0), (2,0), (1,1), (3,1)]:
            if stim == 'antenna':
                return 'radius'
            else:
                return 'width'
        else:
            if stim == 'antenna':
                return 'angle'
            else:
                return 'height'
    else:
        # 2D
        if counter in [0, 2]:
            if stim == 'antenna':
                return 'pos'
            else:
                return 'shape'
        else:
            if stim == 'antenna':
                return 'neg'
            else:
                return 'size'

def dim_label(cond, counter):
    if cond==0:
        if counter in [0, 2]:
            return 'radius'
        else:
            return 'angle'
    elif cond==1:
        if counter in [0, 2]:
            return 'width'
        else:
            return 'height'
    # ii conditions
    elif cond==2:
        if counter in [0, 2]:
            return 'positive'
        else:
            return 'negative'

    elif cond==3:
        if counter in [0, 2]:
            return 'positive (shape)'
        else:
            return 'negative (size)'


def accuracy(sid):
    testtrials = filter(lambda d: d[0]=='test' and d[2]=='feedback' and d[3]=='ncorrect', data(sid))
    return [r[4]/float(N_TEST_TRIALS_PER_BLOCK) for r in testtrials]

def selections(sid, block):
    return np.array([r[5:] for r in filter(lambda d: d[0]=='training' and d[1]==block and d[3]=='selection' and d[4]=='coords', data(sid))])

def selection_feedback(sid, block):
    return [r[-1] for r in filter(lambda d: d[0]=='training' and d[1]==block and d[3]=='feedback', data(sid))]

def plot_selections(sid, blocks=range(8)):
    fig, axi = plt.subplots(1, 8, figsize=(14, 2), sharey=True)
    for block in blocks:
        sel = selections(sid, block)
        lab = np.array(selection_feedback(sid, block))
        chng = np.abs(selection_change(sid, block)).sum(axis=1) >= SEL_CRITERION


        sel_A = sel[np.where((lab=='A') & chng)]
        sel_B = sel[np.where((lab=='B') & chng)]

        rnd_A = sel[np.where((lab=='A') & (~chng))]
        rnd_B = sel[np.where((lab=='B') & (~chng))]


        ax = axi[block]
        if len(sel_A) > 0:
            ax.plot(sel_A[:,0], sel_A[:,1], 's',
                    markeredgecolor='blue', markeredgewidth=1., markerfacecolor='blue', markersize=5)
        if len(sel_B) > 0:
            ax.plot(sel_B[:,0], sel_B[:,1], 's',
                    markeredgecolor='red', markeredgewidth=1., markerfacecolor='red', markersize=5)

        if len(rnd_A) > 0:
            ax.plot(rnd_A[:,0], rnd_A[:,1], 's',
                    markeredgecolor='blue', markeredgewidth=1., markerfacecolor='None', markersize=5, alpha=.3)
        if len(rnd_B) > 0:
            ax.plot(rnd_B[:,0], rnd_B[:,1], 's',
                    markeredgecolor='red', markeredgewidth=1., markerfacecolor='None', markersize=5, alpha=.3)

        ax.set_xlim(-.1, 1.1)
        ax.set_ylim(-.1, 1.1)
    plt.tight_layout()
    plt.show()




def sample_distance(x, rule_type, rule_variant, rescale=False):
    x = np.array(x)

    if rule_type == 'rb':
        if rule_variant in [0, 2]:
            dist = np.abs(0.5 - x[0])
        elif rule_variant in [1, 3]:
            dist = np.abs(0.5 - x[1])

        if rescale:
            dist = dist/MED_SAMPLE_DIST_1D

    elif rule_type == 'ii':

        u = np.array([1., 1.]) / np.sqrt(2) # unit vector

        if rule_variant in [0, 2]:
            dist = np.sqrt(np.sum((x - u * np.dot(x, u)) ** 2))

        elif rule_variant in [1, 3]:

            x2 = np.array([x[0], 1 - x[1]]) # flip to make projection easier
            dist = np.sqrt(np.sum((x2 - u * np.dot(x2, u)) ** 2))

        if rescale:
            dist = dist/MED_SAMPLE_DIST_2D

    return dist

def sample_distance_by_block(sid, rescale=False):
    rule_type = rule_condition(sid)
    rule_variant = rule_counter(sid)
    dist = np.array([[sample_distance(x, rule_type, rule_variant, rescale=rescale) \
                      for x in selections(sid, block)] for block in range(8)])
    return dist

def testdata(sid, block):
    testitems = np.array([r[4:] for r in filter(lambda d: d[0]=='test' and d[1]==block and d[3]=='init_coords', data(sid))], dtype=float)
    testresp = [r[-1] for r in filter(lambda d: d[0]=='test' and d[1]==block and d[3]=='classify' and d[4]=='response', data(sid))]
    testcorrect = [r[-1] for r in filter(lambda d: d[0]=='test' and d[1]==block and d[3]=='classify' and d[4]=='correct', data(sid))]
    return [[float(i[0][0]), float(i[0][1]), i[1], 1*i[2]] for i in zip(testitems, testresp, testcorrect)]

def plot_test_responses(sid, blocks=range(8)):
    fig, axi = plt.subplots(1, 8, figsize=(14, 2), sharey=True)
    for block in blocks:
        td = testdata(sid, block)

        td_A = np.array([x[:2] for x in filter(lambda d: d[2]=='A', td)])
        td_B = np.array([x[:2] for x in filter(lambda d: d[2]=='B', td)])

        ax = axi[block]
        if len(td_A) > 0:
            ax.plot(td_A[:,0], td_A[:,1], 'o', color='blue', markersize=5)
        if len(td_B) > 0:
            ax.plot(td_B[:,0], td_B[:,1], 'o', color='red', markersize=5)
        ax.set_xlim(-.1, 1.1)
        ax.set_ylim(-.1, 1.1)
    plt.tight_layout()
    plt.show()


def selection_change(sid, block):
    init_coord = np.array([r[4:] for r in filter(lambda d: d[0]=='training' and d[1]==block and d[3]=='init_coords', data(sid))])
    sel = selections(sid, block)

    valid = map(lambda r: list(r).count(None)==0, sel)
    d = []
    for i in range(len(sel)):
        if valid[i]:
            d.append(np.array(sel[i], dtype=float) - init_coord[i])
        else:
            d.append(np.array([0., 0.]))

    return np.round(d, 3)


def n_selection_changes(sid):
    return np.transpose([map(lambda f: np.sum(np.abs(f) >= SEL_CRITERION), selection_change(sid, block).transpose()) for block in range(8)])


def n_selection_changes_either_dim(sid):
    return [np.sum(np.sum(np.abs(selection_change(sid, b)), axis=1) > 0.) \
            for b in range(8)]
