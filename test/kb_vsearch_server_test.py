import unittest
import os
import json
import time
import requests  # added

from os import environ
from ConfigParser import ConfigParser
from requests_toolbelt import MultipartEncoder  # added
from pprint import pprint

from biokbase.workspace.client import Workspace as workspaceService
from biokbase.AbstractHandle.Client import AbstractHandle as HandleService  # added

from kb_vsearch.kb_vsearchImpl import kb_vsearch


class kb_vsearchTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = environ.get('KB_AUTH_TOKEN', None)
        cls.ctx = {'token': token, 'provenance': [{'service': 'kb_vsearch',
            'method': 'please_never_use_it_in_production', 'method_params': []}],
            'authenticated': 1}
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('kb_vsearch'):
            cls.cfg[nameval[0]] = nameval[1]
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = workspaceService(cls.wsURL, token=token)
        cls.serviceImpl = kb_vsearch(cls.cfg)

        cls.shockURL = cls.cfg['shock-url']  # added
        cls.handleURL = cls.cfg['handle-service-url']  # added


    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_kb_vsearch_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx


    # Helper script borrowed from the transform service, logger removed
    def upload_file_to_shock(self,
                             shock_service_url = None,
                             filePath = None,
                             ssl_verify = True,
                             token = None):
        """
        Use HTTP multi-part POST to save a file to a SHOCK instance.
        """

        if token is None:
            raise Exception("Authentication token required!")

        #build the header
        header = dict()
        header["Authorization"] = "Oauth {0}".format(token)

        if filePath is None:
            raise Exception("No file given for upload to SHOCK!")

        dataFile = open(os.path.abspath(filePath), 'rb')
        m = MultipartEncoder(fields={'upload': (os.path.split(filePath)[-1], dataFile)})
        header['Content-Type'] = m.content_type

        #logger.info("Sending {0} to {1}".format(filePath,shock_service_url))
        try:
            response = requests.post(shock_service_url + "/node", headers=header, data=m, allow_redirects=True, verify=ssl_verify)
            dataFile.close()
        except:
            dataFile.close()
            raise

        if not response.ok:
            response.raise_for_status()

        result = response.json()

        if result['error']:
            raise Exception(result['error'][0])
        else:
            return result["data"]


    #BEGIN kb_vsearch

    # call this method to get the WS object info of a Paired End Library (will
    # upload the example data if this is the first time the method is called during tests)
    def getSingleEndLibInfo(self,name):
        if hasattr(self.__class__, 'SingleEndLibInfo'):
            return self.__class__.SingleEndLibInfo
        # 1) upload files to shock
        token = self.ctx['token']
        forward_shock_file = self.upload_file_to_shock(
            shock_service_url = self.shockURL,
            filePath = 'data/'+name,
            token = token
            )
        #pprint(forward_shock_file)

        # 2) create handle
        hs = HandleService(url=self.handleURL, token=token)
        forward_handle = hs.persist_handle({
                                        'id' : forward_shock_file['id'], 
                                        'type' : 'shock',
                                        'url' : self.shockURL,
                                        'file_name': forward_shock_file['file']['name'],
                                        'remote_md5': forward_shock_file['file']['checksum']['md5']})

        
        # 3) save to WS
        single_end_library = {
            'lib': {
                'file': {
                    'hid':forward_handle,
                    'file_name': forward_shock_file['file']['name'],
                    'id': forward_shock_file['id'],
                    'url': self.shockURL,
                    'type':'shock',
                    'remote_md5':forward_shock_file['file']['checksum']['md5']
                },
                'encoding':'UTF8',
                'type':'fasta',
                'size':forward_shock_file['file']['size']
            },
            'sequencing_tech':'artificial reads'
        }

        new_obj_info = self.ws.save_objects({
                        'workspace':self.getWsName(),
                        'objects':[
                            {
                                'type':'KBaseFile.SingleEndLibrary',
                                'data':single_end_library,
                                'name':name,
                                'meta':{},
                                'provenance':[
                                    {
                                        'service':'kb_vsearch',
                                        'method':'test_kb_vsearch'
                                    }
                                ]
                            }]
                        })
        self.__class__.SingleEndLibInfo = new_obj_info[0]
        return new_obj_info[0]

    
    def test_VSearch_BasicSearch(self):

        # figure out where the test data lives
        se_lib_info_one = self.getSingleEndLibInfo('input_one.fna')
        pprint(se_lib_info_one)
        se_lib_info_one = self.getSingleEndLibInfo('input_many.fna')
        pprint(se_lib_info_many)

        # Object Info Contents
        # 0 - obj_id objid
        # 1 - obj_name name
        # 2 - type_string type
        # 3 - timestamp save_date
        # 4 - int version
        # 5 - username saved_by
        # 6 - ws_id wsid
        # 7 - ws_name workspace
        # 8 - string chsum
        # 9 - int size
        # 10 - usermeta meta


        # run VSearch_BasicSearch
        workspace_id = se_lib_info_one[7]
        params = {
            'workspace_id': workspace_id,
            'input_one_name': se_lib_info_one[1],
            'input_many_name':  se_lib_info_many[1],
            'output_filtered_name': 'output_filtered.SingleEndLibrary',
            'maxaccepts': 1000,
            'maxrejects': 100000000,
            'wordlength': 8,
            'minwordmatches': 10,
            'ident_thresh': 0.97,
            'ident_mode': 2
        }


        result = self.getImpl().VSearch_BasicSearch(self.getContext(),params)
        print('RESULT:')
        pprint(result)

        # check the output
        info_list = self.ws.get_object_info([{'ref':workspace_id + '/output_filtered.SingleEndLibrary'}], 1)
        self.assertEqual(len(info_list),1)
        output_filtered_info = info_list[0]
        self.assertEqual(output_filtered_info[1],'output_filtered.SingleEndLibrary')
        self.assertEqual(output_filtered_info[2].split('-')[0],'KBaseAssembly.SingleEndLibrary')
        #self.assertEqual(output_filtered_info[10]['Number contigs'],'2')



    #END kb_vsearch

    
