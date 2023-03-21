import os
import ipywidgets as widgets
from IPython.display import display
from pprint import pprint
import ROOT


def build_custom_hypo(hypo_cpp_path,
                      rootlogon_path='/home/hep/packages/kfcmd_tr_ph_v9/share/kfcmd/rootlogon.C'):
    out = widgets.Output(layout={'border': '1px solid black'})
    display(out)
    with out:
        print('Loading macro: {}'.format(rootlogon_path))

    ROOT.gInterpreter.ProcessLine('gROOT->LoadMacro("{}")'.format(rootlogon_path))

    with out:
        print('Building: .L {}++'.format(hypo_cpp_path))

    ROOT.gInterpreter.ProcessLine('.L {}++'.format(hypo_cpp_path))


def run_kinfit(
        input_path,
        output_path='test.root',
        mfield=1.3,
        tr_ph_path='TrPh.C',
        rootlogon_path='/home/hep/packages/kfcmd_tr_ph_v9/share/kfcmd/rootlogon.C'):
    varibles = locals()
    out = widgets.Output(layout={'border': '1px solid black'})
    display(out)
    with out:
        pprint(varibles)

    with out:
        print('Loading macro: {}'.format(rootlogon_path))

    ROOT.gInterpreter.ProcessLine('gROOT->LoadMacro("{}")'.format(rootlogon_path))
    with out:
        print('Building: .L {}++'.format(tr_ph_path))

    ROOT.gInterpreter.ProcessLine(".L {}++".format(tr_ph_path))
    with out:
        print('Opening input file: {}'.format(input_path))

    ROOT.gInterpreter.ProcessLine('TFile fl("{}", "read")'.format(input_path))
    with out:
        print('Creating TrPh object...')

    ROOT.gInterpreter.ProcessLine('TrPh a(tr_ph)')

    with out:
        print('Starting kinfit loop...')

    ROOT.gInterpreter.ProcessLine('a.Loop("{output_path}", {mfield})'.format(output_path=output_path, mfield=mfield))

    with out:
        print('Closing input file...')

    ROOT.gInterpreter.ProcessLine('fl.Close()')
