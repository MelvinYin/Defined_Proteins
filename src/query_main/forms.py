from django import forms
from django.forms.widgets import HiddenInput

class MyForm(forms.Form):
    descr = forms.CharField(label='descr', max_length=100, required=False)
    descr_file = forms.FileField(label='descr_file', required=False)
    pdb_id = forms.CharField(label='pdb_id', max_length=100)
    cid = forms.CharField(label='cid', max_length=100)
    sno = forms.CharField(label='sno', max_length=100)
    aligned_file = forms.FileField(label='aligned_file', required=False)
    matrix_file = forms.FileField(label='matrix_file', required=False)

class ResultsSelect(forms.Form):
    userid = forms.CharField(label='userid', max_length=100,
                             widget=forms.HiddenInput())
    selection = forms.CharField(label='selection', max_length=100, widget =
    forms.HiddenInput())
