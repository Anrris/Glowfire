
class get_pybind_include(object):
    """Helper class to determine the pybind11 include path
    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import imp, os
        try:
            imp.find_module('pybind11')
        except ImportError:
            os.system('pip install pybind11')

        import pybind11
        return pybind11.get_include(self.user)

class get_include(object):
    """Helper class to determine the pybind11 include path
    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """
    def __init__(self,
            from_url,
            to_local,
            extracted_folder,
            extract_method,
            extracted_filename=None
        ):

        self.url = from_url
        self.local = to_local
        self.extracted_filename = extracted_filename
        self.extracted_folder = extracted_folder
        self.extract_method = extract_method

    def __str__(self):
        import pybind11,os,tarfile,zipfile
        from urllib.request import Request, urlopen

        cpp_external_root = self.local
        if not os.path.exists(cpp_external_root):
            os.mkdir(cpp_external_root)

        if self.extracted_filename is None:
            self.extracted_filename = self.url.split('/')[-1]
        filepath = cpp_external_root+self.extracted_filename
        print("Filepath: ", filepath)

        import os
        if not os.path.exists(filepath):
            req = Request(self.url, headers={'User-Agent': 'Mozilla/5.0'})
            tarobj = urlopen(req).read()
            with open(filepath, 'wb') as f:
                f.write(tarobj)


        folderpath = cpp_external_root+self.extracted_folder
        if not os.path.exists(folderpath):
            print('Extracting to:', folderpath)
            if self.extract_method == 'gz':
                tf = tarfile.open(filepath)
                tf.extractall(path=cpp_external_root)
            elif self.extract_method == 'zip':
                with zipfile.ZipFile(filepath, 'r') as zip_ref:
                    zip_ref.extractall(cpp_external_root)

        return folderpath
