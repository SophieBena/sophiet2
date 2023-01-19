from thermoengine import databank
from elasticsearch_dsl import Search,Q

import numpy as np




def test_get_database_client():
    client = databank.DatabaseClient()
    assert client.client is not None

def _test_legacy_search_databse():
    client = databank.DatabaseClient()
    s = Search(using=client.client, index="lepr_experiment_composition").filter(
        Q("match", type='TraceDs') &
        Q("nested", path="entities", query=Q("match", **{'entities.phaseName': 'Clinopyroxene'})) &
        Q("nested", path="entities", query=Q("match", **{'entities.phaseName': 'Liquid'})))

    contains_cpx = []

    for src in s.scan():
        phases = [entity.phaseName for entity in src.entities]

        contains_cpx.append('Clinopyroxene' in phases)

    assert np.all(contains_cpx)
