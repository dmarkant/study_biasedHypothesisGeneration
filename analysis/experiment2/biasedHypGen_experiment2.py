import json
import numpy as np

N_TEST_TRIALS_PER_BLOCK = 32

# expected distance from random selection strategy
RAND_DISTANCE_1D = .25
RAND_DISTANCE_2D = .235

MED_SAMPLE_DIST_1D = .25
MED_SAMPLE_DIST_2D = .217

DATADIR = '../../data_raw/experiment2'

DATA = {}
COND_LABELS = ['1D-REL', '1D-ABS', '2D-REL', '2D-ABS']



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
    c = filter(lambda d: d[0]=='rule_cond', data(sid))[0][1]
    if c == 'rb':
        return '1D'
    else:
        return '2D'

def dim_mapping(sid):
    return filter(lambda d: d[0]=='dim_mapping', data(sid))[0][1]

def stim_condition(sid):
    return filter(lambda d: d[0]=='stim_cond', data(sid))[0][1]

def rule_counter(sid):
    return filter(lambda d: d[0]=='counterbalance', data(sid))[0][1]

def rule_offset(sid):
    return filter(lambda d: d[0]=='rule_offset', data(sid))[0][1]

def rule_direction(sid):
    c = rule_counter(sid)

    if rule_condition(sid)=='1D':
        if c in [0, 3]:
            return 'negative'
        else:
            return 'positive'
    else:
        if c in [0, 2]:
            return 'comp'
        elif c==1:
            return 'positive'
        else:
            return 'negative'

def rule_relation(sid):
    c = rule_counter(sid)
    if rule_condition(sid)=='1D':
        m = dim_mapping(sid)
        if c in [0, 2]:
            if m==0:
                return '1D-C'
            else:
                return '1D-F'
        else:
            if m==0:
                return '1D-F'
            else:
                return '1D-C'
    else:
        if c in [0, 2]:
            return '2D-P'
        else:
            return '2D-N'

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
    else:
        if counter in [0, 2]:
            return 'pos'
        else:
            return 'neg'

def accuracy(sid):
    return [r[4]/32. for r in filter(lambda d: d[0]=='test' and d[2]=='feedback' and d[3]=='ncorrect', data(sid))]

def selections(sid, block):
    return np.array([r[5:] for r in filter(lambda d: d[0]=='training' and d[1]==block and d[3]=='selection' and d[4]=='coords', data(sid))])

def selections_feature_values(sid, block):
    return np.array([r[5:] for r in filter(lambda d: d[0]=='training' and d[1]==block and d[3]=='selection' and d[4]=='fvalue', data(sid))])

def selection_feedback(sid, block):
    return [r[-1] for r in filter(lambda d: d[0]=='training' and d[1]==block and d[3]=='feedback', data(sid))]

def plot_selections(sid, blocks=range(8), featurescale=False):
    fig, axi = plt.subplots(1, 8, figsize=(14, 2), sharey=True)
    for block in blocks:
        if featurescale:
            sel = selections_feature_values(sid, block)
        else:
            sel = selections(sid, block)
        lab = np.array(selection_feedback(sid, block))
        chng = np.abs(selection_change(sid, block)).sum(axis=1) > 0

        sel_A = sel[np.where((lab=='failure') & chng)]
        sel_B = sel[np.where((lab=='success') & chng)]

        rnd_A = sel[np.where((lab=='failure') & (~chng))]
        rnd_B = sel[np.where((lab=='success') & (~chng))]

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

        if not featurescale:
            ax.set_xlim(-.1, 1.1)
            ax.set_ylim(-.1, 1.1)
    plt.tight_layout()
    plt.show()

OFFSET_1D = .1
OFFSET_2D = .1056

def classify(pt, ruletype, rulecounter, ruleoffset):

    x, y = pt

    if ruletype=='1D':
        if rulecounter == 0:
            if (x < (.5 + OFFSET_1D * ruleoffset)):
                return 1
        if rulecounter == 1:
            if (y > (.5 + OFFSET_1D * ruleoffset)):
                return 1
        if rulecounter == 2:
            if (x > (.5 + OFFSET_1D * ruleoffset)):
                return 1
        if rulecounter == 3:
            if (y < (.5 + OFFSET_1D * ruleoffset)):
                return 1
    elif ruletype=='2D':
        if rulecounter == 0:
            if (y > (x + OFFSET_2D * ruleoffset)):
                return 1
        if rulecounter == 1:
            if ((x + OFFSET_2D * ruleoffset) > (1 - y)):
                return 1
        if rulecounter == 2:
            if ((x + OFFSET_2D * ruleoffset) > y):
                return 1
        if rulecounter == 3:
            if ((x + OFFSET_2D * ruleoffset) < (1-y)):
                return 1

    return 0


def centerpoint(ruletype, rulecounter, ruleoffset):

    if ruletype=='1D':
        if rulecounter == 0:
            x = .5 + OFFSET_1D * ruleoffset
            y = .5
        if rulecounter == 1:
            x = .5
            y = .5 + OFFSET_1D * ruleoffset
        if rulecounter == 2:
            x = .5 + OFFSET_1D * ruleoffset
            y = .5
        if rulecounter == 3:
            x = .5
            y = .5 + OFFSET_1D * ruleoffset
    elif ruletype=='2D':
        d = OFFSET_2D * ruleoffset
        if rulecounter == 0:
            x = .5 - (d/2.)
            y = .5 + (d/2.)
        if rulecounter == 1:
            x = .5 - (d/2.)
            y = .5 - (d/2.)
        if rulecounter == 2:
            x = .5 - (d/2.)
            y = .5 + (d/2.)
        if rulecounter == 3:
            x = .5 - (d/2.)
            y = .5 - (d/2.)
    return np.array([x, y])

def sample_distance(x, ruletype, rulecounter, ruleoffset, rescale=True):
    x = np.array(x)

    if ruletype == '1D':
        b = .5 + OFFSET_1D * ruleoffset
        if rulecounter in [0, 2]:
            dist = np.abs(b - x[0])
        elif rulecounter in [1, 3]:
            dist = np.abs(b - x[1])

        if rescale:
            dist = dist/MED_SAMPLE_DIST_1D

    elif ruletype == '2D':


        u = np.array([1., 1.]) / np.sqrt(2) # unit vector

        if rulecounter in [0, 2]:
            x_off = np.array([x[0] + OFFSET_2D * ruleoffset, x[1]])
            dist = np.sqrt(np.sum((x_off - u * np.dot(x_off, u)) ** 2))

        elif rulecounter in [1, 3]:

            x2 = np.array([x[0], 1 - x[1]]) # flip to make projection easier
            x_off = np.array([x2[0] + OFFSET_2D * ruleoffset, x2[1]])
            dist = np.sqrt(np.sum((x2 - u * np.dot(x2, u)) ** 2))

        if rescale:
            dist = dist/MED_SAMPLE_DIST_2D
    return dist


def sample_distance_by_block(sid):
    rule_type = rule_condition(sid)
    rule_variant = rule_counter(sid)
    rule_off = rule_offset(sid)
    dist = np.array([[sample_distance(x, rule_type, rule_variant, rule_off) \
                      for x in selections(sid, block)] for block in range(8)])
    return dist

def testdata(sid, block, featurescale=False):

    if featurescale:
        testitems = np.array([r[4:] for r in filter(lambda d: d[0]=='test' and d[1]==block and d[3]=='init_fvalue', data(sid))], dtype=float)
    else:
        testitems = np.array([r[4:] for r in filter(lambda d: d[0]=='test' and d[1]==block and d[3]=='init_coords', data(sid))], dtype=float)

    testresp = [r[-1] for r in filter(lambda d: d[0]=='test' and d[1]==block and d[3]=='classify' and d[4]=='response', data(sid))]
    testcorrect = [r[-1] for r in filter(lambda d: d[0]=='test' and d[1]==block and d[3]=='classify' and d[4]=='correct', data(sid))]
    return [[float(i[0][0]), float(i[0][1]), i[1], 1*i[2]] for i in zip(testitems, testresp, testcorrect)]

def plot_test_responses(sid, blocks=range(8), featurescale=False):
    fig, axi = plt.subplots(1, 8, figsize=(14, 2), sharey=True)
    for block in blocks:
        td = testdata(sid, block, featurescale=featurescale)

        td_A = np.array([x[:2] for x in filter(lambda d: d[2]=='failure', td)])
        td_B = np.array([x[:2] for x in filter(lambda d: d[2]=='success', td)])

        ax = axi[block]
        ax.plot(td_A[:,0], td_A[:,1], '.', color='blue', markersize=5)
        ax.plot(td_B[:,0], td_B[:,1], 'o', color='red', markersize=5)

        if not featurescale:
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
        return np.transpose([map(lambda f: np.sum(np.abs(f) > 0), selection_change(sid, block).transpose()) for block in range(8)])

def n_selection_changes_either_dim(sid):
    return [np.sum(np.sum(np.abs(selection_change(sid, b)), axis=1) > 0.) \
            for b in range(8)]

