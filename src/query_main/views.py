import os
import pickle
import time

from bokeh.embed import components, autoload_static
from django.core.files.storage import FileSystemStorage
from django.shortcuts import render
from django.views.decorators.cache import never_cache

from config import paths
from descr import descr_main
from matchers.structure.s_matcher import Matcher as s_Matcher
from matchers.descriptor.d_matcher import Matcher as d_Matcher
from preprocess import find_motif_positions
from query_main.forms import MyForm, ResultsSelect
from site_main import settings
from ui.ui_manager import MatcherManager, SegmentCalculator

hash_map = dict()

descr_map = dict()
descr_map['ef'] = "EF-Hand"
descr_map['gxgxxg'] = "GxGxxG"
descr_map['gxxgxg'] = "GxxGxG"
descr_map['custom'] = "Custom"

@never_cache
def results1(request):
    userid = request.POST['userid']
    selection = int(request.POST['selection'])
    global hash_map
    pdb_id, cid, sno, descr, s1, d1, s2, d2, s3, d3, s4, d4, s5, d5 = hash_map[userid]
    if selection == 1:
        html_tag = 'results1.html'
    elif selection == 2:
        html_tag = 'results2.html'
    elif selection == 3:
        html_tag = 'results3.html'
    else:
        form = MyForm()
        return render(request, 'query_page.html', {'form': form})
    return render(request, html_tag, {'pdb_id': pdb_id, 'cid': cid, 'sno': sno,
                                      'descr': descr_map[descr],
                                      'userid': userid, 'd1': d1, 's1': s1,
                                      'd2': d2, 's2': s2, 'd3': d3, 's3': s3,
                                      'd4': d4, "s4": s4, 'd5': d5, 's5': s5})

@never_cache
def results2(request):
    userid = request.POST['userid']
    selection = int(request.POST['selection'])
    global hash_map
    pdb_id, cid, sno, descr, s1, d1, s2, d2, s3, d3, s4, d4, s5, d5 = \
        hash_map[userid]
    if selection == 1:
        html_tag = 'results1.html'
    elif selection == 2:
        html_tag = 'results2.html'
    elif selection == 3:
        html_tag = 'results3.html'
    else:
        form = MyForm()
        return render(request, 'query_page.html', {'form': form})
    return render(request, html_tag, {'pdb_id': pdb_id, 'cid': cid, 'sno': sno,
                                      'descr': descr_map[descr],
                                      'userid': userid, 'd1': d1, 's1': s1,
                                      'd2': d2, 's2': s2, 'd3': d3, 's3': s3,
                                      'd4': d4, "s4": s4, 'd5': d5, 's5': s5})

@never_cache
def results3(request):
    userid = request.POST['userid']
    selection = int(request.POST['selection'])
    global hash_map
    pdb_id, cid, sno, descr, s1, d1, s2, d2, s3, d3, s4, d4, s5, d5 = \
        hash_map[userid]
    if selection == 1:
        html_tag = 'results1.html'
    elif selection == 2:
        html_tag = 'results2.html'
    elif selection == 3:
        html_tag = 'results3.html'
    else:
        form = MyForm()
        return render(request, 'query_page.html', {'form': form})
    return render(request, html_tag, {'pdb_id': pdb_id, 'cid': cid, 'sno': sno,
                                      'descr': descr_map[descr],
                                      'userid': userid, 'd1': d1, 's1': s1,
                                      'd2': d2, 's2': s2, 'd3': d3, 's3': s3,
                                      'd4': d4, "s4": s4, 'd5': d5, 's5': s5})

@never_cache
def index(request):
    if request.method == 'POST':
        form = MyForm(request.POST, request.FILES)
        media_dir = os.path.join(paths.INTERNAL, "media")
        if request.FILES and ('aligned_file' in request.FILES or
                              "matrix_file" in request.FILES):
            pdb_seq_file = paths.RCSB_SEQS_FASTA
            motif_file = os.path.join(media_dir, "motif.txt")
            if 'aligned_file' in request.FILES:
                aligned_file = request.FILES['aligned_file']
                fs = FileSystemStorage(location=media_dir)
                filename = fs.save(aligned_file.name, aligned_file)
                aligned_file = os.path.join(media_dir, filename)
                assert os.path.isfile(pdb_seq_file)
                assert os.path.isfile(aligned_file)
                find_motif_positions.from_aligned(pdb_seq_file, aligned_file,
                                                  pdb_seq_file, motif_file,
                                                  key='custom')
            else:
                matrix_file = request.FILES['matrix_file']
                fs = FileSystemStorage(location=media_dir)
                filename = fs.save(matrix_file.name, matrix_file)
                matrix_file = os.path.join(media_dir, filename)
                assert os.path.isfile(matrix_file)
                assert os.path.isfile(pdb_seq_file)
                find_motif_positions.from_nbdb(matrix_file, motif_file,
                                               pdb_seq_file, 0, key="custom")
            assert os.path.isfile(motif_file)
            with open(motif_file, 'rb') as file:
                motif_pos_map = pickle.load(file)
            descrs = descr_main.calculate(motif_pos_map)

            s_matcher = s_Matcher(cropped=False)
            s_matcher.load(descrs)
            d_matcher = d_Matcher()
            d_matcher.load(descrs)

            s_map = dict()
            s_map['custom'] = s_matcher
            d_map = dict()
            d_map['custom'] = d_matcher

            segment_calculater = SegmentCalculator(matchers=s_map)

            descriptor_file_html = os.path.join("static/tmp",
                                                "descriptor.pkl")
            descriptor_file_self = os.path.join(settings.BASE_DIR,
                                                descriptor_file_html)

            with open(descriptor_file_self, 'wb') as file:
                pickle.dump((s_map, d_map, segment_calculater), file, -1)

            return render(request, 'results4.html', {'descr_file':
                                                         descriptor_file_html})

        elif form.is_valid():
            pdb_id = request.POST['pdb_id']
            cid = request.POST['cid']
            sno = request.POST['sno']
            descr = request.POST['descr']

            if not descr and not request.FILES:
                form = MyForm()
                x = render(request, 'query_page.html', {'form': form})
                x["Cache-Control"] = "no-cache, no-store, must-revalidate"  #
                # HTTP 1.1.
                x["Pragma"] = "no-cache"  # HTTP 1.0.
                return x
            if request.FILES:
                descr = 'custom'
                descr_file = request.FILES['descr_file']
                fs = FileSystemStorage(location=media_dir)
                filename = fs.save(descr_file.name, descr_file)
                full_path = os.path.join(media_dir, filename)
                manager = MatcherManager(descr, pdb_id, cid, sno, full_path)
            else:
                manager = MatcherManager(descr, pdb_id, cid, sno)
            userid = str(hash(time.time()))
            global hash_map
            s1, d1 = components(manager.ui.e1)
            s2, d2 = components(manager.ui.e2)
            s3, d3 = components(manager.ui.e3)
            s4, d4 = components(manager.ui.e4)
            s5, d5 = components(manager.ui.e5)

            hash_map[userid] = [pdb_id, cid, sno, descr, s1, d1, s2, d2, s3,
                                d3, s4, d4, s5, d5]
            html_tag = 'results1.html'
            return render(request, html_tag,
                          {'pdb_id': pdb_id, 'cid': cid, 'sno': sno,
                           'descr': descr_map[descr], 'userid': userid,
                           'd1': d1, 's1': s1, 'd2': d2, 's2': s2, 'd3': d3,
                           's3': s3, 'd4': d4, "s4": s4, 'd5': d5, 's5': s5})
    else:
        form = MyForm()
    x = render(request, 'query_page.html', {'form': form})
    x["Cache-Control"] = "no-cache, no-store, must-revalidate"  # HTTP 1.1.
    x["Pragma"] = "no-cache"  # HTTP 1.0.
    return x