from django import forms
from pseudogene.models import *


class NewForm(forms.ModelForm):
    class Meta:
        model = NewModel
        fields = ['field1', 'field2', 'field3']

    def clean_field1(self):
        field1 = self.cleaned_data['field1']

        return field1 + '_add'


class AnotherForm(NewForm):
    class Meta(NewForm.Meta):
        model = NewModel
        exclude = ['field3']


class AllForm(forms.Form):
    field1 = forms.CharField()
    field2 = forms.CharField()
    field3 = forms.CharField()

    def clean(self):
        cleaned_data = super(AllForm, self).clean()
        field1 = cleaned_data.get('field1')

        return field1 + '_add'

    def __init__(self, *args, **kwargs):
        super(AllForm, self).__init__(*args, **kwargs)
        for field_name, field in self.fields.items():
            field.widget.attrs['placeholder'] = 'Fill in...'


class Form1(AllForm):
    field2 = None
    field4 = forms.CharField()




