# import earthchem
from elasticsearch import Elasticsearch
# from elasticsearch_dsl import Search,Q

class DatabaseClient:
    def __init__(self):
        self._init_login_info()
        self._init_client()


    def _init_login_info(self):
        self.username = 'lepr1'
        self.password = 'Earthchemlepr1!'
        self.kibana = 'https://elasticsearch-staging.earthchem.org/_plugin/kibana'
        self.host = 'https://elasticsearch-staging.earthchem.org'

    def _init_client(self, scheme="https", port=443,
                     use_ssl=True, verify_certs=True, timeout=60,
                     max_retries=10, retry_on_timeout=True):

        client = Elasticsearch(
            [self.host],
            http_auth=(self.username, self.password),
            scheme=scheme,
            port=port,
            use_ssl=use_ssl,
            verify_certs=verify_certs,
            timeout=timeout,
            max_retries=max_retries,
            retry_on_timeout=retry_on_timeout
        )
        self.client = client