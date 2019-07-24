import django_tables2 as tables


class NewTable(tables.Table):
    field2 = tables.Column()
    field3 = tables.Column()
    field4 = tables.Column()

    class Meta:
        attrs = {'id': 'table1'}

    @staticmethod
    def make_table(a_list):
        data = [{
            'field2': a_list[0],
            'field3': a_list[1],
            'field4': 4
        }]
        table = NewTable(data)
        return table


class AnotherTable(NewTable):
    field1 = tables.Column()
    field5 = tables.Column(verbose_name=' ')

    class Meta:
        attrs = {'id': 'table1'}
        sequence = ('field1', 'field2')

    @staticmethod
    def make_table(a_list):
        data = [{
            'field1': a_list[0],
            'field2': a_list[1],
            'field3': a_list[2],
            'field4': 4,
            'field5': 'g'
        }]
        table = AnotherTable(data)
        return table

