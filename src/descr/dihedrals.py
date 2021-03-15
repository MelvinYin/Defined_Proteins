from descr import geometry


def get_descr_dihedrals(C, CA, N, dsr_snos):
    """
    C, CA, N should already be sorted by sno.
    """
    angles = dict()
    dihedrals = _pdbTorsion(C, CA, N)
    phis, psis = list(zip(*dihedrals))
    angles['sno'] = dsr_snos
    angles['phi'] = [360]
    angles['psi'] = [360]
    angles['region'] = [""]
    angles['ss'] = [""]
    for i in range(len(dsr_snos)-2):
        phi = phis[i]
        psi = psis[i]
        region, ss = get_ramachandran_region(phi, psi)
        angles['phi'].append(phi)
        angles['psi'].append(psi)
        angles['region'].append(region)
        angles['ss'].append(ss)
    angles['phi'].append(360)
    angles['psi'].append(360)
    angles['region'].append("")
    angles['ss'].append("")

    assert len(angles['psi']) == len(angles['sno'])
    assert len(angles['phi']) == len(angles['sno'])
    assert len(angles['region']) == len(angles['sno'])
    assert len(angles['ss']) == len(angles['sno'])
    return angles, CA

def _pdbTorsion(C, CA, N):
    dihedrals = []
    for i in range(1, len(C)-1):
        prev_C = C[i-1]
        curr_C = C[i]
        curr_CA = CA[i]
        curr_N = N[i]
        next_N = N[i+1]

        try:
            dihedral = geometry.calcDihedrals(prev_C, curr_N, curr_CA,
                                              curr_C, next_N)
        except ValueError:
            err = f"Dihedral calculation failed."
            raise Exception(err)
        dihedrals.append(dihedral)
    return dihedrals

##############################################################################
# From ramachandran file.
##############################################################################

def get_ramachandran_region(phi, psi):
    if not (-180 <= phi <= 180 and -180 <= psi <= 180):
        raise ValueError("Incorrect value of phi/psi torsion angles")

    phi = int(phi)
    psi = int(psi)

    gap = 10    # grid size in degrees

    # DATA CODE   , 'oBbAaLlexyzw',
    # DATA NEWCOD , 'XXB b A a L l p ~b~a~l~p',
    # DATA RTYPE  ,  1,4,3,4,3,4,3,3,2,2,2,2,

    code = {'o': 1, 'B': 2, 'b': 3, 'A': 4, 'a': 5, 'L': 6, 'l': 7, 'e': 8,
            'x': 9, 'y': 10, 'z': 11, 'w': 12}

    symbol = {'o': 'X', 'B': 'B', 'b': 'b', 'A': 'A', 'a': 'a', 'L': 'L',
              'l': 'l', 'e': 'p', 'x': '~b', 'y': '~a', 'z': '~l', 'w': '~p'}

    rama_map = ['bBBBBBBBBBBbbbxxooooowwwwwoooooxxxxb',
                'bBBBBBBBBBBBbbbxxooooooooooooooxxxbb',
                'bBBBBBBBBBBBbbbxxxxooooooooooooxxbbb',
                'bbBBBBBBBBBBBbbbbxxooooooooooooxxxbb',
                'bbBBBBBBBBBBBBbbxxxooooooooooooxxxxb',
                'bbBBBBBBBBBBBbbbxxxoooooooooooooxxxb',
                'bbbBBBBBBBBBbbbbxxxooozzzzzzoooooxxb',
                'bbbbBBBBBBBbbbbbbxxzzzzzzzzzoooooxxb',
                'bbbbbBbbBbbbbbbbbxxzzzzzllzzoooooxxb',
                'bbbbbbbbbbbbbbxxxxxzzllllzzzoooooxxx',
                'bbbbbbbbbbbbxxxxxxxzzllllzzzoooooxxx',
                'xbbbbbbbbbbbbxxoooozzzlllzzzzooooxxx',
                'xxbbbbbbbbbbbxxoooozzllllllzzooooxxx',
                'yyaaaaaaaaaayyyoooozzllLllzzzooooxxx',
                'yaaaaaaaaaaayyyoooozzzlLLlzzzooooxxx',
                'yaaaaaaAaaaaayyyooozzzlllllzzooooxxx',
                'yaaaaaAAAAaaayyyyooozzzlllzzzooooxxx',
                'yaaaaaAAAAAaaayyyooozzzzlllzzooooxxx',
                'yaaaaAAAAAAAaaayyyoozzzlzllzzooooxxx',
                'yaaaaaAAAAAAAaayyyyozzzzzzzzzooooxxx',
                'yyaaaaAAAAAAAAaayyyozzzzzzzzzooooxxx',
                'yyaaaaaAAAAAAAaaayyyoooooooooooooxxx',
                'yyyaaaaaAAAAAAAaayyyoooooooooooooxxx',
                'yaaaaaaaaAAAAAAaaayyyooooooooooooxxx',
                'aayaaaaaaaaAAAAaaayyyooooooooooooxxx',
                'yyyaaaaaaaaaaaaaaaayyooooooooooooxxx',
                'yyyyyaaaaaaaaaaaaayyyooooooooooooxxx',
                'oyyyyaaaaaaaaaayyyyyyooooooooooooooo',
                'oooyyyyyyyyayyyyyyyyoooooooooooooooo',
                'oooxxxxbbxxxxxxxxooooooooooooooooooo',
                'xxxxxbbxxxxxxxooooooowwwwwoooooooooo',
                'xxxxxxbbbbxxxxooooooowwwwwoooooooooo',
                'xbbbxbbbbbxxxxxoooooowwewwwooooooooo',
                'xbbbbbbbbbbxxxxoooooowwewwwooooooxxx',
                'xbbbbbbbbbbbbxxoooooowweewwooooooxxx',
                'bbbbbbbbbbbbbxxoooooowwwwwwooooooxxb']

    x = (phi + 180) // gap
    x = min(x, 35)
    y = (psi + 180) // gap
    y = min(y, 35)
    m = rama_map[35-y][x]
    return code[m], symbol[m]
