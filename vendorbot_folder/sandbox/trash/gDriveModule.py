import logging, os, pathlib

from pydrive.drive import GoogleDrive
from pydrive.auth import GoogleAuth

class PathTester:
    def __init__(self):
        pass
    def pathman(self):
        # print(os.listdir())
        # print(os.path.abspath(__file__))
        # print(pathlib.Path(__file__).parent.absolute())
        # print(os.path.dirname(__file__))
        # print(os.path.dirname(os.path.abspath(__file__)))
        # print(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "/secrets", "/client_secrets.json")
        # print(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "secrets", "client_secrets.json"))
        # path = pathlib.Path(os.path.abspath(__file__)).parent
        # print(path.parent)
        # print(os.getcwd())
        path = pathlib.Path(os.path.abspath(__file__))
        print(os.path.join((path.parent.parent), "secrets", "client_secrets.json"))

class GoogleDriver:
    def __init__(self, *args, **kwargs):
        if args: # If args is not empty.
            self.args = args
        if kwargs:
            for key, value in kwargs.items():
                setattr(self, key, value)

    def connect(self):
        path = pathlib.Path(os.path.abspath(__file__)).parent.parent # /vendorbot_container/ directory
        GoogleAuth.DEFAULT_SETTINGS['client_config_file'] = os.path.join(path, "secrets", "client_secrets.json") 
        self.gauth = GoogleAuth()
        
        # Try to load saved client credentials
        self.gauth.LoadCredentialsFile(os.path.join(path, "secrets", "mycreds.txt"))

        if self.gauth.credentials is None:
            # Authenticate if they're not there
            # This is what solved the issues:
            self.gauth.GetFlow()
            self.gauth.flow.params.update({'access_type': 'offline'})
            self.gauth.flow.params.update({'approval_prompt': 'force'})
            self.gauth.LocalWebserverAuth()

        elif self.gauth.access_token_expired:
            # Refresh them if expired
            self.gauth.Refresh()
        else:
            # Initialize the saved creds
            self.gauth.Authorize()

        # Save the current credentials to a file
        self.gauth.SaveCredentialsFile(os.path.join(path, "secrets", "mycreds.txt"))
        return True 
    
    def upload (self, name, filepath):
        self.drive = GoogleDrive(self.gauth)
        file1 = self.drive.CreateFile({'title': name})  
        # find a way how to add current data to the file name
        file1.SetContentFile(filepath)
        file1.Upload()
        logging.info("file uploaded")
        return True
        