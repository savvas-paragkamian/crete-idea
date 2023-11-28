#!/usr/bin/env python3
import clarivate.wos_starter.client
from clarivate.wos_starter.client.apis.tags import documents_api
from clarivate.wos_starter.client.model.documents_list import DocumentsList
from pprint import pprint
# Defining the host is optional and defaults to http://api.clarivate.com/apis/wos-starter
# See configuration.py for a list of all supported configuration parameters.
configuration = clarivate.wos_starter.client.Configuration(
    host = "http://api.clarivate.com/apis/wos-starter"
)

# Enter a context with an instance of the API client
with clarivate.wos_starter.client.ApiClient(configuration) as api_client:
    # Create an instance of the API class
    api_instance = documents_api.DocumentsApi(api_client)

    # example passing only required values which don't have defaults set
    query_params = {
        'q': "PY=2020",
    }
    try:
        # Query Web of Science documents 
        api_response = api_instance.documents_get(
            query_params=query_params,
        )
        pprint(api_response)
    except clarivate.wos_starter.client.ApiException as e:
        print("Exception when calling DocumentsApi->documents_get: %s\n" % e)

    # example passing only optional values
    query_params = {
        'db': "WOS",
        'q': "PY=2020",
        'limit': 10,
        'page': 1,
        'sortField': "sortField_example",
    }
    try:
        # Query Web of Science documents 
        api_response = api_instance.documents_get(
            query_params=query_params,
        )
        pprint(api_response)
    except clarivate.wos_starter.client.ApiException as e:
        print("Exception when calling DocumentsApi->documents_get: %s\n" % e)
