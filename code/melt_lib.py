import RNA, numpy as np


def native_pairs(seq, salt, Tref=25):
    md = RNA.md(); md.temperature = Tref; md.salt = salt
    s, g = RNA.fold_compound(seq, md).mfe()
    pt = RNA.ptable(s)
    return s, [(i, pt[i]) for i in range(1, len(seq) + 1) if pt[i] > i]


def fnat(seq, native, T, salt, delta):
    md = RNA.md(); md.temperature = T; md.salt = salt
    fc = RNA.fold_compound(seq, md); n = len(seq)
    if delta:
        for i in range(1, n + 1):
            for j in range(i + 4, n + 1):
                fc.sc_add_bp(i, j, delta)
    s, g = fc.mfe(); fc.exp_params_rescale(g); fc.pf()
    b = fc.bpp()
    return float(np.mean([b[i][j] for i, j in native]))


def tm(Ts, F):
    for k in range(1, len(Ts)):
        if F[k - 1] >= 0.5 > F[k]:
            return Ts[k - 1] + (F[k - 1] - 0.5) * (Ts[k] - Ts[k - 1]) / (F[k - 1] - F[k])
    return None
