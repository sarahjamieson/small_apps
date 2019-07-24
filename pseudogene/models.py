from __future__ import unicode_literals

from django.db import models


class NewModel(models.Model):
    field1 = models.CharField(max_length=10)
    field2 = models.CharField(max_length=10)
    field3 = models.CharField(max_length=10, null=True)

    def compare(self, obj):
        diffs = {}
        excluded_keys = ['_state', 'id']
        d1, d2 = self.__dict__, obj.__dict__
        for k, v in d1.items():
            if k in excluded_keys:
                continue
            if v and d2[k] and v != d2[k]:
                diffs[k] = (str(v), str(d2[k]))

        return diffs
